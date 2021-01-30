# PrimeSieve
Super fast, dynamically expanding prime sieve for fast primality queries, forward or backward iteration.

Queries automatically grow the sieve as needed.

Uses approximately `n/16` bytes of storage, where `n` is the largest value queried/grown to, so for example, a sieve grown to encompass all 32-bit integers (2^32) uses just 256 MB of storage (2^32 / 16 bytes).

Computing up to 2^32 like this typically takes less than half a second.

### Examples of usage: ###

- - - -

Create a sieve but do not initialize.

`PrimeSieve sieve;`

- - - -

Query primality:

`std::cout << (sieve.IsPrime(127) ? "prime" : "composite") << std::endl;`

> prime

- - - -

Get next prime (and not including):

`std::cout << sieve.NextPrime(127) << std::endl;`

> 131

- - - -

Get previous prime (and not including):

`std::cout << sieve.PrevPrime(127) << std::endl;`

> 113

- - - -

Iterate forever upward through primes:

```
for (const U64 p : sieve)
{
    std::cout << p << ", ";
}
```

> 2, 3, 5, 7, 11, 13, 17, 19, 23, 29...

- - - -

Iterate forever upward through primes, starting after (and not including):

```
for (const U64 p : sieve.IterateForwardFrom(50))
{
    std::cout << p << ", ";
}
```

> 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101...

- - - -

Iterate downward through primes, starting before (and not including):

```
for (const U64 p : sieve.IterateBackwardFrom(50))
{
    std::cout << p << ", ";
}
```

> 47, 43, 41, 37, 31, 29, 23, 19, 17, 13, 11, 7, 5, 3, 2,

- - - -

Initialize prime sieve by pre-computing up to a value `x` using efficient multi-threaded evaluation, such that subsequent queries <= `x` are constant-time lookups and don't require growing the sieve:

`PrimeSieve sieve(x);`

- - - -

Grow sieve efficiently up to value `x`:

`sieve.GrowTo(x);`

- - - -
