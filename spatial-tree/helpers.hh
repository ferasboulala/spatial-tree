#pragma once

#include <tuple>

#include "compilation.hh"

namespace st {

template <typename SequenceIterator>
struct EnumerationIterator {
    SequenceIterator it;
    uint64_t         idx;

    __always_inline EnumerationIterator(SequenceIterator beg) : it(beg), idx(0) {}

    __always_inline EnumerationIterator& operator++() {
        ++it;
        ++idx;

        return *this;
    }

    __always_inline std::tuple<uint64_t, decltype(*it)> operator*() const {
        return std::make_tuple(idx, *it);
    }

    __always_inline bool operator==(const EnumerationIterator<SequenceIterator>& other) {
        return it == other.it;
    }
    __always_inline bool operator!=(const EnumerationIterator<SequenceIterator>& other) {
        return it != other.it;
    }
};

template <typename SequenceT>
struct EnumerationContainer {
    const SequenceT& sequence;

    using sequence_iterator_t = decltype(sequence.begin());

    __always_inline EnumerationContainer(const SequenceT& seq) : sequence(seq) {}

    __always_inline EnumerationIterator<sequence_iterator_t> begin() const {
        return EnumerationIterator<sequence_iterator_t>(sequence.begin());
    }

    __always_inline EnumerationIterator<sequence_iterator_t> end() const {
        return EnumerationIterator(sequence.end());
    }
};

template <typename SequenceT>
__always_inline auto enumerate(const SequenceT& sequence) {
    return EnumerationContainer(sequence);
}

}  // namespace st