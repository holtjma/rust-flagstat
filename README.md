# rust-flagstat
Place for building a tool similar to [samtools flagstat](http://www.htslib.org/doc/samtools-flagstat.html)
for the purpose of understanding how the numbers are generated and to try out [rust-htslib](https://crates.io/crates/rust-htslib).
This testing tool does not perfectly replicate the output format.
However, the numbers are the same for reads where all reads are marked as "passed".

Note: there possibly some edge cases that are not captured by my test data that will create slightly different statistics.
For example, the test data used did not contain any duplicates, so they were not present in our test set for validating the relevant numbers.
