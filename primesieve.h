/*******************************************************************
*
*    Author: Kareem Omar
*    kareem.h.omar@gmail.com
*    https://github.com/komrad36
*
*    Last updated Jan 29, 2021
*******************************************************************/

#pragma once

#include <cstdint>
#include <immintrin.h>

class PrimeSieve
{
    // using decls
    using U32 = uint32_t;
    using U64 = uint64_t;

    // forward decls
private:
    class FwdIterator;
    class FwdIteratorEndSentinel;
    class FwdIteratorFrom;
    class RevIterator;
    class RevIteratorEndSentinel;
    class RevIteratorFrom;

    // constants
private:
    static constexpr U64 kBitsPerSeg = 3ULL * 5 * 7 * 11 * 13 * 17;
    static constexpr U64 kBlocksPerSeg = (kBitsPerSeg + 63) >> 6;
    static constexpr U64 kUnusedBitsPerSeg = 64 - (kBitsPerSeg & 63);
    static constexpr U64 kMaxThreads = 32;

    // static methods
private:
    static U64 ComputeAutoNumThreads();
    static void ComputeInternal(U64* __restrict const pSieve, U64 iStart, U64 iEnd);

    static U64 CountLeadingZeros(U64 x)
    {
        return _lzcnt_u64(x);
    }

    static U64 CountTrailingZeros(U64 x)
    {
        return _tzcnt_u64(x);
    }

    static U64 ShiftRight(U64 x, U32 i)
    {
#if defined(__clang__) || defined(__GNUC__)
        return x >> (i & 63U);
#else
        return _shrx_u64(x, i);
#endif
    }

    static U64 ShiftRight(U64 x, U64 i)
    {
#if defined(__clang__) || defined(__GNUC__)
        return x >> (i & 63ULL);
#else
        return _shrx_u64(x, static_cast<U32>(i));
#endif
    }

    static U64 ShiftLeft(U64 x, U32 i)
    {
#if defined(__clang__) || defined(__GNUC__)
        return x << (i & 63U);
#else
        return _shlx_u64(x, i);
#endif
    }

    static U64 ShiftLeft(U64 x, U64 i)
    {
#if defined(__clang__) || defined(__GNUC__)
        return x << (i & 63ULL);
#else
        return _shlx_u64(x, static_cast<U32>(i));
#endif
    }

    static U64 SmartRem(U64 a, U64 b)
    {
        // clang has a bug that results in suboptimal
        // codegen if you try this, but it generates
        // the right thing anyway if you don't
#ifdef __clang__
        return a % b;
#else
        return static_cast<U32>((a | b) >> 32) ? a % b : static_cast<U32>(a) % static_cast<U32>(b);
#endif
    }

    // public methods
public:

    // x: initial value to compute sieve up to, i.e. so that primality queries <= x will then be constant time
    // 0 == no initial computation (default)
    //
    // numThreads: thread count for large computations
    // 0 == auto (default)
    PrimeSieve(U64 x = 0, U64 numThreads = 0) : m_numThreads(numThreads ? numThreads > kMaxThreads ? kMaxThreads : numThreads : ComputeAutoNumThreads())
    {
        GrowTo(x);
    }

    ~PrimeSieve();

    // compute sieve up to x, i.e. so that primality queries for <= x will then be constant time
    void GrowTo(U64 x)
    {
        if (x >= 3)
            GrowToInternal((x - 1) / (kBitsPerSeg << 1) + 1);
    }

    // query primality of x, automatically growing sieve as necessary to reach x
    bool IsPrime(U64 x)
    {
        if (x == 2)
            return true;

        if (!(x & 1))
            return false;

        const U64 iSeg = x / (kBitsPerSeg << 1);
        if (iSeg >= m_numSegsComputed)
            GrowToInternal(iSeg + 1);

        x = (x >> 1) + kUnusedBitsPerSeg * iSeg;

        return !(m_pSieve[x >> 6] & (1ULL << (x & 63)));
    }

    // returns an object which iterates up through primes,
    // starting at the next prime after (and not including) x
    // and continues forever
    FwdIteratorFrom IterateForwardFrom(U64 x)
    {
        return FwdIteratorFrom(*this, x);
    }

    // returns an object which iterates down through primes,
    // starting at the previous prime before (and not including) x
    // and continues until all primes are visited (i.e. ends after visiting 2)
    RevIteratorFrom IterateBackwardFrom(U64 x)
    {
        return RevIteratorFrom(*this, x);
    }

    // return next prime after (and not including) x
    U64 NextPrime(U64 x)
    {
        return *IterateForwardFrom(x).begin();
    }

    // return the previous prime before (and not including) x
    // behavior is undefined for x <= 2
    U64 PrevPrime(U64 x)
    {
        return *IterateBackwardFrom(x).begin();
    }

    FwdIterator begin()
    {
        return *this;
    }

    FwdIteratorEndSentinel end()
    {
        return FwdIteratorEndSentinel();
    }

    // forward iterator
private:
    class FwdIteratorEndSentinel {};

    class FwdIterator
    {
        friend class PrimeSieve;

        FwdIterator(PrimeSieve& sieve) : m_sieve(sieve), m_pSieve(m_sieve.m_pSieve), m_block(0), m_iBlock(~0ULL), m_iEndBlock(PrimeSieve::kBlocksPerSeg * m_sieve.m_numSegsComputed), m_x(2) {}
        FwdIterator(PrimeSieve& sieve, U64 x) : m_sieve(sieve)
        {
            const U64 iBit_NativeSpace = (x + 1) >> 1;
            const U64 iSeg = iBit_NativeSpace / PrimeSieve::kBitsPerSeg;
            if (iSeg >= sieve.m_numSegsComputed)
                m_sieve.GrowToInternal(iSeg + 1);

            const U64 iBit_SegSpace = iBit_NativeSpace + PrimeSieve::kUnusedBitsPerSeg * iSeg;
            m_iEndBlock = PrimeSieve::kBlocksPerSeg * m_sieve.m_numSegsComputed;
            m_pSieve = m_sieve.m_pSieve;
            m_iBlock = iBit_SegSpace >> 6;
            m_block = (~m_pSieve[m_iBlock]) & (~0ULL << (iBit_SegSpace & 63));
            AdvanceInternal();
        }

        void AdvanceInternal()
        {
            while (!m_block)
            {
                if (++m_iBlock >= m_iEndBlock)
                {
                    m_sieve.GrowToInternal(m_sieve.m_numSegsComputed + 1);
                    m_iEndBlock = PrimeSieve::kBlocksPerSeg * m_sieve.m_numSegsComputed;
                    m_pSieve = m_sieve.m_pSieve;
                }

                m_block = ~m_pSieve[m_iBlock];
            }

            m_x = (m_iBlock << 7) + (CountTrailingZeros(m_block) << 1) - (kUnusedBitsPerSeg << 1) * (m_iBlock / kBlocksPerSeg) + 1;
        }

    public:
        bool operator==(const PrimeSieve::FwdIteratorEndSentinel&)
        {
            return false;
        }

        bool operator!=(const PrimeSieve::FwdIteratorEndSentinel&)
        {
            return true;
        }

        void operator++()
        {
            m_block &= m_block - 1;
            AdvanceInternal();
        }

        U64 operator*() const
        {
            return m_x;
        }

    private:
        PrimeSieve& m_sieve;
        const U64* __restrict m_pSieve;
        U64 m_block;
        U64 m_iBlock;
        U64 m_iEndBlock;
        U64 m_x;
    };

    class FwdIteratorFrom
    {
        friend class PrimeSieve;

        FwdIteratorFrom(PrimeSieve& sieve, U64 x) : m_sieve(sieve), m_x(x) {}

    public:
        FwdIterator begin()
        {
            return m_x >= 2 ? FwdIterator(m_sieve, m_x) : FwdIterator(m_sieve);
        }

        FwdIteratorEndSentinel end()
        {
            return FwdIteratorEndSentinel();
        }

    private:
        PrimeSieve& m_sieve;
        U64 m_x;
    };

    // reverse iterator
private:
    class RevIteratorEndSentinel {};

    class RevIterator
    {
        friend class PrimeSieve;

        RevIterator(PrimeSieve& sieve, U64 x) : m_sieve(sieve)
        {
            if (x <= 2)
            {
                m_x = 1;
                return;
            }

            const U64 iBit_NativeSpace = x >> 1;
            const U64 iSeg = iBit_NativeSpace / PrimeSieve::kBitsPerSeg;
            if (iSeg >= sieve.m_numSegsComputed)
                m_sieve.GrowToInternal(iSeg + 1);

            const U64 iBit_SegSpace = iBit_NativeSpace + PrimeSieve::kUnusedBitsPerSeg * iSeg;
            m_pSieve = m_sieve.m_pSieve;
            m_iBlock = iBit_SegSpace >> 6;
            m_block = ~(m_pSieve[m_iBlock] | ShiftLeft(~0ULL, iBit_SegSpace));
            m_x = x;
            AdvanceInternal();
        }

        void AdvanceInternal()
        {
            while (!m_block)
            {
                if (m_iBlock == 0)
                {
                    --m_x;
                    return;
                }
                m_block = ~m_pSieve[--m_iBlock];
            }

            m_x = (m_iBlock << 7) + ((CountLeadingZeros(m_block) ^ 63) << 1) - (kUnusedBitsPerSeg << 1) * (m_iBlock / kBlocksPerSeg) + 1;
        }

    public:
        bool operator==(const PrimeSieve::RevIteratorEndSentinel&)
        {
            return m_x == 1;
        }

        bool operator!=(const PrimeSieve::RevIteratorEndSentinel&)
        {
            return m_x != 1;
        }

        void operator++()
        {
            m_block &= ~ShiftRight(1ULL << 63, CountLeadingZeros(m_block));
            AdvanceInternal();
        }

        U64 operator*() const
        {
            return m_x;
        }

    private:
        PrimeSieve& m_sieve;
        const U64* __restrict m_pSieve;
        U64 m_block;
        U64 m_iBlock;
        U64 m_x;
    };

    class RevIteratorFrom
    {
        friend class PrimeSieve;

        RevIteratorFrom(PrimeSieve& sieve, U64 x) : m_sieve(sieve), m_x(x) {}

    public:
        RevIterator begin()
        {
            return RevIterator(m_sieve, m_x);
        }

        RevIteratorEndSentinel end()
        {
            return RevIteratorEndSentinel();
        }

    private:
        PrimeSieve& m_sieve;
        U64 m_x;
    };

    // private methods
private:
    void GrowToInternal(U64 newNumSegs);
    void ComputeToInternal(U64 newNumSegs);

    // member variables
private:
    U64* __restrict m_pSieve = nullptr;
    U64 m_numSegsAllocated = 0;
    U64 m_numSegsComputed = 0;
    U64 m_numThreads;
};
