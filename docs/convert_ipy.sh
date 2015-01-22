 #!/bin/sh

 find . -name '*.ipynb' -execdir ipython nbconvert --to rst {} \;§
