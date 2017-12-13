# Package install

## Inital setup
- Clone this Github repository:
```bash
git clone https://github.com/seignovert/shdompp.git
```

- Add the bin folder to your `$PATH`:
```bash
cd shdompp
echo "# SHDOM bin folder:" >> $HOME/.bashrc
echo "export PATH=$(pwd)/bin:\$PATH" >> $HOME/.bashrc
source $HOME/.bashrc
```

- Compile SHDOMPP scripts:
```bash
make
```

## Testing
Testing scripts (with there input data) are available in the `tests` directory.

All the tests can be run with:
```bash
make test
```

_Note:_ For now `run_compare_shdompp.sh` does not work out of the box.
