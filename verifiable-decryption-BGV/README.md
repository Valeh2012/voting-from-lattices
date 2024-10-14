# verifiable-decryption-BGV

The code is forked from [tjesi/verifiable-decryption-BGV](https://github.com/tjesi/verifiable-decryption-BGV) repository. Parameters are adjusted in `parameters.h` header file to support BGV encryption with ~62 bit modulus. Make sure [NTL](https://libntl.org/doc/tour-unix.html) is installed in your system.

First, uncomment correct parameter sets in `parameters.h`, and then run `make run` to see results. For your own parameters, use `params.py` to calculate required standard deviations and bound values.