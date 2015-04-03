 #!/bin/sh

 find . -name '*.ipynb' -execdir ipython nbconvert --to rst {} \;

find . -name '*.ipynb' -execdir ipython nbconvert --to python {} \;
