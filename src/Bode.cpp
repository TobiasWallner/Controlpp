#include <controlpp/Bode.hpp>

namespace controlpp{

    std::ostream& operator<<(std::ostream& stream, EBodeCsvReadError val){
        const char* str = "";
        switch(val){
            case EBodeCsvReadError::CouldNotFindFrequencyVector: str = "Could not find frequency vector"; break;
            case EBodeCsvReadError::CouldNotFindAmplitudeVectors: str = "Could not find amplitude vector"; break;
        };
        return stream << str;
    }


    namespace{
        /**
         * \brief makes an ASCII character lowercase
         */
        char to_lowercase(char c){
            switch(c){
                case 'A': return 'a';
                case 'B': return 'b';
                case 'C': return 'c';
                case 'D': return 'd';
                case 'E': return 'e';
                case 'F': return 'f';
                case 'G': return 'g';
                case 'H': return 'h';
                case 'I': return 'i';
                case 'J': return 'j';
                case 'K': return 'k';
                case 'L': return 'l';
                case 'M': return 'm';
                case 'N': return 'n';
                case 'O': return 'o';
                case 'P': return 'p';
                case 'Q': return 'q';
                case 'R': return 'r';
                case 'S': return 's';
                case 'T': return 't';
                case 'U': return 'u';
                case 'V': return 'v';
                case 'W': return 'w';
                case 'X': return 'x';
                case 'Y': return 'y';
                case 'Z': return 'z';
                default: return c;
            }
        }
        
        /**
         * \brief makes a string lowercase in place
         */
        void to_lowercase_ip(std::string& str){
            for(auto& elem : str){
                elem = to_lowercase(elem);
            }
        }
    }

    std::expected<Bode<double>, std::variant<EBodeCsvReadError, csvd::ReadError>> read_bode_from_csv(
        std::istream& stream, 
        const csvd::Settings& csv_settings,
        EFrequencyInterpretation freq_interp,
        EMagnitudeInterpretation mag_interp,
        EPhaseInterpretation phase_interp
    ){
        std::expected<csvd::CSVd, csvd::ReadError> csv_result = csvd::read(stream, csv_settings);
        if(csv_result.has_value()){

            csvd::CSVd& csv = csv_result.value();
            // good case

            // make all names lowercase
            for(csvd::Column& column : csv){
                to_lowercase_ip(column.name);
            }

            // find frequency
            const auto freq_column = csv.find_if([](std::string_view name){return name.starts_with("f");});
            if(freq_column == csv.end()){
                return std::unexpected(EBodeCsvReadError::CouldNotFindFrequencyVector);
            }

            // determine frequency unit
            bool freq_hz = false;
            switch(freq_interp){
                case EFrequencyInterpretation::AutoHz:{
                    const bool found_rad = freq_column->name.contains("rad");
                    if(found_rad){
                        freq_hz = false;
                    }else{
                        freq_hz = true;
                    }
                }break;
                case EFrequencyInterpretation::AutoRad:{
                    const bool found_hz = freq_column->name.contains("hz");
                    if(found_hz){
                        freq_hz = true;
                    }else{
                        freq_hz = false;
                    }
                }break;
                case EFrequencyInterpretation::ForceHz:{
                    freq_hz = true;
                }break;
                case EFrequencyInterpretation::ForceRad:{
                    freq_hz = false;
                }break;
            }

            // create frequency vector
            Eigen::VectorXd frequs;
            frequs.resize(freq_column->data.size());
            if(freq_hz){
                for(size_t i = 0; i != freq_column->data.size(); ++i){
                    frequs[i] = freq_column->data[i];
                }
            }else{
                for(size_t i = 0; i != freq_column->data.size(); ++i){
                    frequs[i] = freq_column->data[i] / (2 * std::numbers::pi);
                }
            }
            

            // find real and imaginary
            const auto real_column = csv.find_if([](std::string_view name){return name.starts_with("re");});
            const auto imag_column = csv.find_if([](std::string_view name){return name.starts_with("im");});

            if(real_column != csv.end() && imag_column != csv.end()){
                // found real and imaginary
                Eigen::VectorXcd values;
                values.resize(real_column->data.size());
                for(size_t i = 0; i != real_column->data.size(); ++i){
                    values[i].real(real_column->data[i]);
                    values[i].imag(imag_column->data[i]);
                }
                Bode bode(std::move(frequs), std::move(values));
                return bode;
            }

            // find magnitude and phase
            const auto mag_column = csv.find_if([](std::string_view name){return name.starts_with("mag");});
            const auto phase_column = csv.find_if([](std::string_view name){return name.starts_with("ph");});
        
            if(mag_column != csv.end() && phase_column != csv.end()){
                // found magnitude and phase

                // determine magnitude unit

                // determine magnitude unit
                bool is_mag_dB = false;
                switch(mag_interp){
                    case EMagnitudeInterpretation::Auto : {
                        is_mag_dB = mag_column->name.contains("db");
                    }break;
                    case EMagnitudeInterpretation::ForceAbs : {
                        is_mag_dB = false;
                    }break;
                    case EMagnitudeInterpretation::ForceDB : {
                        is_mag_dB = true;
                    }break;
                };

                // create magnitude vector with absolute values
                Eigen::VectorXd mags;
                mags.resize(mag_column->data.size());
                if(is_mag_dB){
                    for(size_t i = 0; i != mag_column->data.size(); ++i){
                        mags[i] = std::pow(10.0, mag_column->data[i] / 20.0);
                    }
                }else{
                    for(size_t i = 0; i != mag_column->data.size(); ++i){
                        mags[i] = mag_column->data[i];
                    }
                }
                
                
                // determine phase unit
                bool phase_is_rad = true;
                switch(phase_interp){
                    case EPhaseInterpretation::AutoRad : {
                        phase_is_rad = true;
                        const bool found_deg = phase_column->name.contains("deg");
                        if(found_deg){
                            phase_is_rad = false;
                        }
                    }break;
                    case EPhaseInterpretation::AutoDeg : {
                        phase_is_rad = false;
                        const bool found_rad = phase_column->name.contains("rad");
                        if(found_rad){
                            phase_is_rad = true;
                        }
                    }break;
                    case EPhaseInterpretation::ForceRad : {
                        phase_is_rad = true;
                    }break;
                    case EPhaseInterpretation::ForceDeg : {
                        phase_is_rad = false;
                    }break;
                }

                // make phase vector in rad
                Eigen::VectorXd phases;
                phases.resize(phase_column->data.size());
                if(phase_is_rad){
                    for(size_t i = 0; i != phase_column->data.size(); ++i){
                        phases[i] = phase_column->data[i];
                    }
                }else{
                    for(size_t i = 0; i != phase_column->data.size(); ++i){
                        phases[i] = phase_column->data[i] * (2.0 * std::numbers::pi);
                    }
                }

                // create complex value vector
                Eigen::VectorXcd values;
                {
                    auto real = mags.array() * phases.array().cos();
                    auto imag = mags.array() * phases.array().sin();
                    values.real() = real;
                    values.imag() = imag;
                }

                Bode bode(std::move(frequs), std::move(values));
                return bode;
            }
            return std::unexpected(EBodeCsvReadError::CouldNotFindAmplitudeVectors);
        }else{
            // error case
            return std::unexpected(csv_result.error());
        }
    }

}// namespace controlpp