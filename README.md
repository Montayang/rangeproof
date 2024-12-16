# rangeproof

It is the final version(maybe, with zk) of the protocol of polynomial rangeproof.

It is based on mcl, used the MSM function provided by mcl

Some tools are implemented:

- uni-veriate NTT and INTT, bi-veriate NTT and INTT(no-recursion version)
- uni-veriate polynomial product, bi-veriate polynomial product
- uni-veriate polynomial division, bi-veriate polynomial division(not optimal)
- uni-veriate kzg, bi-veriate kzg and batching version of each of them
- sumset representation

Some baseline protocol implemented:
- Plookup
- Bulletproof

Each protocol outputs the prover time, the verifier time and proof size.
