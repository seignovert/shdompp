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
```bash
cd script
./run_shdompp_examples
```
