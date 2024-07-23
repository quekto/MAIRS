# Description

Reads all `.spa`/`.SPA` files within current directory or deeper and converts them into `.csv`.

The `read_spa()` function is from [lerkoah / spa-on-python](https://github.com/lerkoah/spa-on-python.git), which is based on [LoadSprecta](https://la.mathworks.com/matlabcentral/fileexchange/57904-loadspectra) function from matlab.

## Usage

When reading all files within and below the current directory:

```terminal
python convert_csv_to_spa.py
```

Target directory can be specified:

```terminal
python convert_csv_to_spa.py target/directory/
```

## License

GNU General Public License v3.0 or later

See [LICENSE](https://github.com/quekto/MAIRS/blob/main/LICENSE) to see the full text.
