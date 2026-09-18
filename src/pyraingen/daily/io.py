"""Reading the generator's inputs.

Replaces the Fortran ``readdata`` and the ``nearby_station_details.out``
reader.  ``data_r.dat`` is gone: parameters are passed as arguments now.
``nearby_station_details.out`` is still read from disk because
:func:`pyraingen.getnearbystations.station` writes it and users edit it.
"""

import os

import numpy as np

from .calendarutil import NDAY, NOUT, daycount

MAX_YEARS = 200


def read_nearby_stations(path, nstn):
    """Read the first ``nstn`` rows of a ``nearby_station_details.out`` file.

    Returns
    -------
    tuple of ndarray
        ``(idx, weight)`` -- station numbers and their raw (un-cumulated)
        weights.

    Notes
    -----
    The file also carries year counts, start years and mean annual rainfall,
    but the Fortran overwrote all three from the daily data files, so they are
    ignored here too.
    """
    with open(path) as fh:
        lines = fh.readlines()

    rows = lines[6:]
    idx = np.zeros(nstn + 1, dtype=np.int64)
    weight = np.zeros(nstn + 1, dtype=np.float32)

    found = 0
    for line in rows:
        parts = line.split()
        if len(parts) < 6:
            continue
        found += 1
        idx[found] = int(parts[1])
        weight[found] = np.float32(float(parts[2]))
        if found == nstn:
            break

    if found < nstn:
        raise ValueError(
            f"{path} lists {found} nearby stations, but nstation={nstn} was requested"
        )
    return idx, weight


def read_station_data(data_path, station_ids, nstn, nout=NOUT):
    """Read the ``rev_drNNNNN.txt`` daily files into an ``rf`` array.

    Mirrors the Fortran ``readdata``: the record starts at the first 1 January
    encountered, a trailing partial year is dropped, and the mean annual
    rainfall is the record total divided by the number of whole years.

    Returns
    -------
    tuple
        ``(rf, nyrs, nstrt, avrf)`` with ``rf`` shaped
        ``(nstn+1, nyear_max+1, 13, 32)`` in float32 and the rest 1-based.
    """
    per_station = []
    nyrs = np.zeros(nstn + 1, dtype=np.int64)
    nstrt = np.zeros(nstn + 1, dtype=np.int64)
    avrf = np.zeros(nstn + 1, dtype=np.float32)

    for jj in range(1, nstn + 1):
        # The Fortran built this name by overwriting digits 8-12 of the
        # literal 'rev_dr046115.txt', so the leading zero is fixed and the
        # station index is always five digits.
        fname = os.path.join(data_path, f"rev_dr0{int(station_ids[jj]):05d}.txt")
        years, months, days, amounts = _read_daily_file(fname)
        per_station.append((years, months, days, amounts))

    nyear_max = 1
    parsed = []
    for jj in range(1, nstn + 1):
        years, months, days, amounts = per_station[jj - 1]
        start, count, rows = _lay_out_record(years, months, days, nout)
        nyrs[jj] = count
        nstrt[jj] = start
        nyear_max = max(nyear_max, count)
        parsed.append(rows)

    if nyear_max > MAX_YEARS:
        raise ValueError(
            f"a station has {nyear_max} years of record, above the {MAX_YEARS}-year limit"
        )

    rf = np.zeros((nstn + 1, nyear_max + 2, 13, 32), dtype=np.float32)
    for jj in range(1, nstn + 1):
        years, months, days, amounts = per_station[jj - 1]
        slot = parsed[jj - 1]
        keep = slot >= 1
        rf[jj, slot[keep], months[keep], days[keep]] = amounts[keep]

    for jj in range(1, nstn + 1):
        avrf[jj] = _mean_annual(rf[jj], int(nyrs[jj]), int(nstrt[jj]), nout)

    return rf, nyrs, nstrt, avrf


def _read_daily_file(fname):
    """Parse one ``year month day amount`` file, skipping its header line.

    The Fortran read used ``err=107``, which abandoned the rest of the file on
    the first unparseable line rather than skipping it.  That is reproduced
    here so a corrupt file truncates the record the same way.
    """
    years = []
    months = []
    days = []
    amounts = []
    with open(fname) as fh:
        fh.readline()
        for line in fh:
            parts = line.split()
            if len(parts) < 4:
                if not parts:
                    continue
                break
            try:
                iy = int(parts[0])
                im = int(parts[1])
                idy = int(parts[2])
                amt = float(parts[3])
            except ValueError:
                break
            years.append(iy)
            months.append(im)
            days.append(idy)
            amounts.append(amt)

    if not years:
        raise ValueError(f"{fname} contains no usable daily records")

    return (
        np.array(years, dtype=np.int64),
        np.array(months, dtype=np.int64),
        np.array(days, dtype=np.int64),
        np.array(amounts, dtype=np.float32),
    )


def _lay_out_record(years, months, days, nout):
    """Assign each row a 1-based year slot, reproducing ``readdata``'s rules.

    The record opens at the first 1 January in the file; rows before it are
    discarded.  A new slot starts whenever the calendar year increases.  If the
    file does not end on 31 December, the final partial year is dropped.
    """
    slot = np.zeros(years.shape, dtype=np.int64)

    first = np.flatnonzero((months == 1) & (days == 1))
    if first.size == 0:
        raise ValueError("no 1 January found: cannot establish the record start")
    begin = first[0]
    nstart = int(years[begin])

    count = 0
    current = nstart - 1
    for n in range(begin, years.size):
        iy = int(years[n])
        if iy < nstart:
            continue
        if iy > current:
            count += 1
            current = iy
        slot[n] = count

    # readdata: `if(i2.ne.nout.and.i3.ne.31) i1=i1-1` -- a trailing partial
    # year is dropped only when the last row is neither in December nor on a
    # 31st.  Reproduced exactly, quirk included.
    if months[-1] != nout and days[-1] != 31:
        count -= 1

    slot[slot > count] = 0
    return nstart, count, slot


def _mean_annual(rf_stn, nyear, nstart, nout):
    """Mean annual rainfall over the whole-year part of the record."""
    total = np.float32(0.0)
    for i in range(1, nyear + 1):
        for j in range(1, nout + 1):
            nd = int(NDAY[j])
            if j == 2:
                nd = daycount(nstart - 1, i)
            total = np.float32(total + rf_stn[i, j, 1 : nd + 1].sum(dtype=np.float32))
    if nyear <= 0:
        return np.float32(0.0)
    return np.float32(total / np.float32(nyear))
