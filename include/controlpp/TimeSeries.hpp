#pragma once

#include <ostream>
#include <Eigen/Dense>

namespace controlpp
{
    /**
     * \brief Contiains time and values pairs
     * 
     * Stores two Eigen vectors: `.times` and `.values` that store the time series and are puplicly accessible.
     * 
     * \tparam T The datatype used for the `times` and `values` vectors
     * \tparam N The size of the timeseries. Can be a static size or `Eigen::Dynamic` to determine the size at runtime.
     */
    template<class T, int N = Eigen::Dynamic>
    struct TimeSeries{
        Eigen::Vector<T, N> times;
        Eigen::Vector<T, N> values;
    };

    /**
     * \brief Contiains time and values pairs
     * 
     * Overload that for the case with runtime size using `N=Eigen::Dynamic` as a template parameter (=default).
     * Offers constructors to set the size of the vectors.
     * 
     * Stores two Eigen vectors: `.times` and `.values` that store the time series and are puplicly accessible.
     * 
     * \tparam T The datatype used for the `times` and `values` vectors
     * \tparam N The size of the timeseries. Can be a static size or `Eigen::Dynamic` to determine the size at runtime.
     */
    template<class T>
    class TimeSeries<T, Eigen::Dynamic>{
    public:
        Eigen::Vector<T, Eigen::Dynamic> times;
        Eigen::Vector<T, Eigen::Dynamic> values;

        /// @brief Constructor that sets the size of the timeseries
        /// @param size The size being allocated for the `times` and `values` vectors
        TimeSeries(int size)
            : times(size)
            , values(size){}

        /// @brief Default copy constructor
        TimeSeries(const TimeSeries&) = default;
        /// @brief Default move constructor 
        TimeSeries(TimeSeries&&) = default;
        /// @brief Default copy assignment 
        TimeSeries& operator=(const TimeSeries&) = default;
        /// @brief Default move assignment 
        TimeSeries& operator=(TimeSeries&&) = default;
    };

    /// @brief Prints the timeseries as a `.csv` file
    /// @tparam T The datatype of the timeseries
    /// @tparam N The size of the timeseries
    /// @param stream The output stream/file that should be printed to
    /// @param timeseries The timeseries holding the data
    /// @return a reference to the stream for `operator<<` chaining
    template<class T, int N>
    std::ostream& operator<< (std::ostream& stream, const TimeSeries<T, N>& timeseries){
        stream << "times, values\n";
        for(int i = 0; (i < timeseries.times.size()) && (i < timeseries.values.size()); ++i){
            stream << timeseries.times(i) << ", " << timeseries.values(i) << "\n";
        }
        return stream;
    }
} // namespace controlpp
