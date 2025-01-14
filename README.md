# rangeproof

It is the final version(maybe, with zk) of the protocol of polynomial rangeproof.

It is based on mcl, used the MSM function provided by mcl

Some tools are implemented:

- uni-variate NTT and INTT, bi-variate NTT and INTT(no-recursion version)
- uni-variate polynomial product, bi-variate polynomial product
- uni-variate polynomial division, bi-variate polynomial division(not optimal)
- uni-variate kzg, bi-variate kzg and batching version of each of them
- sumset representation

Some baseline protocol implemented:
- Plookup
- Bulletproof

Each protocol outputs the prover time, the verifier time and proof size.
