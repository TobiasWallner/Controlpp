#pragma once

// std
#include <utility>
#include <optional>

// Eigen
#include <Eigen/Dense>

namespace controlpp{


    /**
     * \brief Finds elements in a range that enclose v.
     * 
     * Searches for two elements in the range given by `[first, last)` such that `*itr <= v && v <= *(itr+1)` is true.
     * 
     * \param first the first iterator of the range (points to the first element of the range)
     * \param last the last iterator of the range (points past the last element of the range)
     * \param v the input value to search for
     * \returns a pair of `[low, high]` iterators that enclose the input value `v` or `std::nullopt` if no enclosing sub-range could be found.
     */
    template<class Itr, class T>
    std::optional<std::pair<Itr, Itr>> find_enclosing(Itr first, Itr last, const T& v){
        Itr itr_low = first;
        Itr itr_high = first;
        ++itr_high;

        for(; itr_high != last; (void)++itr_low, (void)++itr_high){
            if(*itr_low <= v && v <= *itr_high){
                std::pair<Itr, Itr> result(itr_low, itr_high);
                return result;
            }
        }

        return std::nullopt;
    }

    /**
     * \brief Finds elements in a range that enclose v.
     * \param range The range to search in/iterate through
     * \param v The value to search for
     * \returns a pair of `[low, high]` iterators that enclose the input value `v` or `std::nullopt` if no enclosing sub-range could be found.
     * \see template<class Itr, class T> std::optional<std::tuple<Itr, Itr>> find_enclosing(Itr first, Itr last, const T& v)
     */
    template<class T, int Size = Eigen::Dynamic>
    std::optional<std::pair<const T*, const T*>> find_enclosing(const Eigen::Vector<T, Size>& range, const T& v){
        const T* first = range.data();
        const T* last = range.data() + range.size();
        return find_enclosing(first, last, v);
    }

}