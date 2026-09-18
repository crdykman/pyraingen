# pyraingen

A package for stochastically generating daily and subdaily rainfall in Australia with ifd constraining.

## Installation

```bash
$ pip install pyraingen
```

## Usage

`pyraingen` can be used to stochastically generate regionalised daily rainfall, disaggregate daily rainfall to subdaily fragments and constrain generated rainfall to observed or predicted Intensity Frequency Duration (IFD) relationships.
The three main functions are:

```python
from pyraingen.regionaliseddailysim import regionaliseddailysim
from pyraingen.regionalisedsubdailysim import regionalisedsubdailysim
from pyraingen.ifdcond import ifdcond
```

See [`docs/example.ipynb`](docs/example.ipynb) for a worked example of all three.

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`pyraingen` was created by Caleb Dykman. Caleb Dykman retains all rights to the source and it may not be reproduced, distributed, or used to create derivative works.

## Credits

`pyraingen` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).
