#pragma once

#include <ostream>
#include <Eigen/Dense>

namespace controlpp
{

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
    class TimeSeries{
    private:
        Eigen::Vector<T, Eigen::Dynamic> times_;
        Eigen::Vector<T, Eigen::Dynamic> values_;

    public:

        TimeSeries() = default;

        /**
         * @brief Initialises an empty (= uninitialised) timeseries with a given allocation size (number of value-time pairs)
         * 
         * Only allocates the storage. You have to fill the storage yourself
         * 
         * @param size The number of value time pairs that this timeseries holds
         * @see resize(size_t)
         */
        TimeSeries(int size) : times_(size), values_(size){}

        /**
         * @brief 
         * @param t Time vector (Needs to be ordered ascending)
         * @param v Value vector
         */
        TimeSeries(
            const Eigen::Vector<T, Eigen::Dynamic>& t, 
            const Eigen::Vector<T, Eigen::Dynamic>& v)
            : times_(t)
            , values_(v){}

        /**
         * @brief 
         * @param t Time vector (Needs to be ordered ascending)
         * @param v Value vector
         */
        TimeSeries(
            Eigen::Vector<T, Eigen::Dynamic>&& t, 
            Eigen::Vector<T, Eigen::Dynamic>&& v)
            : times_(std::move(t))
            , values_(std::move(v)){}

        TimeSeries(
            const Eigen::Vector<T, Eigen::Dynamic>& t, 
            Eigen::Vector<T, Eigen::Dynamic>&& v)
            : times_(t)
            , values_(std::move(v)){}

        TimeSeries(
            Eigen::Vector<T, Eigen::Dynamic>&& t, 
            const Eigen::Vector<T, Eigen::Dynamic>& v)
            : times_(std::move(t))
            , values_(v){}

        Eigen::Vector<T, Eigen::Dynamic>& times() {return this->times_;}
        const Eigen::Vector<T, Eigen::Dynamic>& times() const {return this->times_;}

        T& times(size_t i){return this->times_(i);}
        const T& times(size_t i) const {return this->times_(i);}

        Eigen::Vector<T, Eigen::Dynamic>& values() {return this->values_;}
        const Eigen::Vector<T, Eigen::Dynamic>& values() const {return this->values_;}

        T& values(size_t i){return this->values_(i);}
        const T& values(size_t i) const {return this->values_(i);}

        /// @brief Default copy constructor
        TimeSeries(const TimeSeries&) = default;
        /// @brief Default move constructor 
        TimeSeries(TimeSeries&&) = default;
        /// @brief Default copy assignment 
        TimeSeries& operator=(const TimeSeries&) = default;
        /// @brief Default move assignment 
        TimeSeries& operator=(TimeSeries&&) = default;

        void resize(int n){
            this->times_.resize(n);
            this->values_.resize(n);
        }

        size_t size() const {return std::min(this->times_.size(), this->values_.size());}
    };

    /// @brief Prints the timeseries as a `.csv` file
    /// @tparam T The datatype of the timeseries
    /// @tparam N The size of the timeseries
    /// @param stream The output stream/file that should be printed to
    /// @param timeseries The timeseries holding the data
    /// @return a reference to the stream for `operator<<` chaining
    template<class T>
    std::ostream& operator<< (std::ostream& stream, const TimeSeries<T>& timeseries){
        stream << "Times (s), Values\n";
        for(size_t i = 0; (i < timeseries.size()) && (i < timeseries.size()); ++i){
            stream << timeseries.times(i) << ", " << timeseries.values(i) << "\n";
        }
        return stream;
    }

    template<class T>
    class ComplexTimeSeries{
    public:
        Eigen::Vector<T, Eigen::Dynamic> times;
        Eigen::Vector<std::complex<T>, Eigen::Dynamic> values;

        ComplexTimeSeries() = default;

        ComplexTimeSeries(int size) : times(size), values(size) {}

        ComplexTimeSeries(
            const Eigen::Vector<T, Eigen::Dynamic>& t, 
            const Eigen::Vector<std::complex<T>, Eigen::Dynamic>& v)
            : times(t)
            , values(v){}

        ComplexTimeSeries(
            Eigen::Vector<T, Eigen::Dynamic>&& t, 
            Eigen::Vector<std::complex<T>, Eigen::Dynamic>&& v)
            : times(std::move(t))
            , values(std::move(v)){}

        ComplexTimeSeries(
            const Eigen::Vector<T, Eigen::Dynamic>& t, 
            Eigen::Vector<std::complex<T>, Eigen::Dynamic>&& v)
            : times(std::move(t))
            , values(std::move(v)){}

        ComplexTimeSeries(
            Eigen::Vector<T, Eigen::Dynamic>&& t, 
            const Eigen::Vector<std::complex<T>, Eigen::Dynamic>& v)
            : times(std::move(t))
            , values(std::move(v)){}

        /// @brief Default copy constructor
        ComplexTimeSeries(const ComplexTimeSeries&) = default;
        /// @brief Default move constructor 
        ComplexTimeSeries(ComplexTimeSeries&&) = default;
        /// @brief Default copy assignment 
        ComplexTimeSeries& operator=(const ComplexTimeSeries&) = default;
        /// @brief Default move assignment 
        ComplexTimeSeries& operator=(ComplexTimeSeries&&) = default;

        void resize(int n){
            times.resize(n);
            values.resize(n);
        }
    };

    template<class T>
    std::ostream& operator<< (std::ostream& stream, const ComplexTimeSeries<T>& timeseries){
        stream << "Times (s), Real, Imag, Abs, Arg\n";
        for(int i = 0; (i < timeseries.times.size()) && (i < timeseries.values.size()); ++i){
            const T t = timeseries.times(i);
            const T real = std::real(timeseries.values(i));
            const T imag = std::imag(timeseries.values(i));
            const T abs = std::abs(timeseries.values(i));
            const T arg = std::arg(timeseries.values(i));
            stream << t << ", " << real << ", " << imag << ", " << abs << ", " << arg << "\n";
        }
        return stream;
    }
    
} // namespace controlpp
