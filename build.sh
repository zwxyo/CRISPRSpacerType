#!/bin/bash

PIP_NO_INDEX='False'

echo "Python version:"
$PYTHON --version

echo "Conda channels:"
conda config --show channels

echo "Pip config:"
pip config list

$PYTHON -m pip install --no-deps --ignore-installed -r $SRC_DIR/config/requirements.txt --no-cache-dir


mkdir -p $PREFIX/share/CRISPRSpacerType
cp -r src config database bin external $PREFIX/share/CRISPRSpacerType

mkdir -p $PREFIX/bin
chmod +x $PREFIX/share/CRISPRSpacerType/bin/*.sh

# ln -s $PREFIX/share/CRISPRSpacerType/bin/CRISPRSpacerType.sh $PREFIX/bin/CRISPRSpacerType
for script in $PREFIX/share/CRISPRSpacerType/bin/*.sh; do
    name=$(basename "$script" .sh)
    ln -s "$script" "$PREFIX/bin/$name"
done