/**
 * @file csv_to_hty.cpp
 * @author Sarutch Supaibulpipat (Pokpong) {ssupaibu@cmkl.ac.th}
 * @brief 
 *  Converts CSV files into HTY binary files with metadata
 * @version 0.1
 * @date 2024-11-29
 */
#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <string>

#include "../third_party/nlohmann/json.hpp"
#include "hty_helper.hpp"

using json = nlohmann::json;

/**
 * @brief Converts a CSV file to an HTY binary file with metadata.
 *
 * @param csv_file_path Path to the input CSV file.
 * @param hty_file_path Path to the output HTY binary file.
 */
void convert_from_csv_to_hty(std::string csv_file_path,
                             std::string hty_file_path)
{
    std::ifstream csv(csv_file_path);
    if (!csv.is_open()) {
        std::cerr << "Cannot open file " << csv_file_path << "\n";
        return;
    }
    std::ofstream hty(hty_file_path, std::ios::binary);
    if (!hty.is_open()) {
        std::cerr << "Cannot open file " << hty_file_path << "\n";
        return;
    }

    try {
        std::string col_info;
        if (!std::getline(csv, col_info)) {
            throw std::runtime_error("empty...");
        }

        StrVec_t col_names = split(col_info, ",");
        trim_all(col_names);
        if (std::any_of(col_names.begin(), col_names.end(),
                        [](std::string &s) { return !has_alpha(s); })) {
            throw std::runtime_error("Invaild Column name/s");
        } else if (!is_unique(col_names)) {
            throw std::runtime_error("Not unique");
        }

        std::list<StrVec_t> rows;

        std::vector<Num32Type> types(col_names.size(), Num32Type::I);

        std::string row;
        while (std::getline(csv, row)) {
            StrVec_t cols = split(row, ",");
            if (cols.size() != col_names.size()) {
                throw std::runtime_error("Column number inconsistent");
            }
            trim_all(cols);
            for (std::size_t i = 0; i < cols.size(); i++) {
                if (str_is_int(cols[i])) {
                } else if (str_is_float(cols[i])) {
                    types[i] = Num32Type::F;
                } else {
                    throw std::runtime_error("Column type inconsistent");
                }
            }
            rows.push_back(cols);
        }

        std::vector<Num32> data(rows.size() * rows.back().size());
        std::vector<Num32>::iterator data_i = data.begin();
        for (StrVec_t cols : rows) {
            for (std::size_t i = 0; i < cols.size(); i++) {
                switch (types[i]) {
                case Num32Type::I:
                    *data_i = std::stoi(cols[i]);
                    break;
                case Num32Type::F:
                    *data_i = std::stof(cols[i]);
                    break;
                }
                data_i++;
            }
        }
        json metadata;
        metadata["num_rows"] = data.size() / col_names.size();
        metadata["num_groups"] = 1;
        metadata["groups"] = json::array();
        metadata["groups"][0] = {{"num_columns", col_names.size()},
                                 {"offset", 0},
                                 {"columns", json::array()}};
        json &col_metadata = metadata["groups"][0]["columns"];
        for (int i = 0; i < col_names.size(); i++) {
            col_metadata[i]["column_name"] = col_names[i];
            col_metadata[i]["column_type"] = num32_type_str_bimap.get(types[i]);
        }
        auto count_char = [](std::uint32_t num) -> std::uint32_t {
            return std::to_string(num).size();
        };
        std::uint32_t offset =
            metadata.dump().size() + sizeof(std::uint32_t) - 1;
        offset += count_char(offset + count_char(offset));
        metadata["groups"][0]["offset"] = offset;
        std::string metadata_str = metadata.dump();
        std::uint32_t metadata_size = metadata_str.size();
        hty.write((char *)&metadata_size, sizeof(metadata_size));
        hty.write(&metadata_str[0], metadata_size);
        hty.write((char *)&data[0], data.size() * sizeof(Num32));
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
    }
    csv.close();
    hty.close();
}

int main()
{
    std::string csv_file_path, hty_file_path;

    // Get arguments
    std::cout << "Please enter the .csv file path:" << std::endl;
    std::cin >> csv_file_path;
    std::cout << "Please enter the .hty file path:" << std::endl;
    std::cin >> hty_file_path;

    // Convert
    convert_from_csv_to_hty(csv_file_path, hty_file_path);

    return 0;
}