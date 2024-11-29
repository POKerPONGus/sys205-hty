/**
 * @file analyze.cpp
 * @author Sarutch Supaibulpipat (Pokpong) {ssupaibu@cmkl.ac.th}
 * @brief
 *
 *  Conatins function for basic querying and adding data to a HTY file.
 * @version 0.1
 * @date 2024-11-29
 */
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "../third_party/nlohmann/json.hpp"
#include "hty_helper.hpp"

typedef std::pair<std::uint32_t, std::uint32_t> UInt32Pair_t;
typedef std::pair<std::uint32_t, UInt32Vec_t> UInt32ValVecPair_t;

struct Num32Vec {
    std::vector<Num32> vec;
    Num32Type type;
};

enum class FilterOp {
    EQ,  // Equal
    L,   // Less than
    G,   // Greater than
    LEQ, // Less than or equal
    GEQ, // Greater than or equal
    NEQ  // Not equal
};

using json = nlohmann::json;

/* Function Delcarations */

template <typename T> static bool (*_get_op(FilterOp op))(T, T);
std::optional<UInt32Pair_t> find_col_idxs(const json &metadata,
                                          std::string target_col);
std::optional<UInt32ValVecPair_t> find_many_col_idx(const json &metadata,
                                                    StrVec_t target_cols);
UInt32Vec_t calc_row_offsets(UInt32Vec_t col_idxs, std::uint32_t total_col);
json update_metadata(json &metadata, std::uint32_t new_row_count);

std::optional<json> extract_metadata(std::string hty_file_path);
std::optional<Num32Vec> project_one(json &metadata, std::string hty_file_path,
                                    std::string projected_column);
std::optional<std::vector<Num32Vec>> project_many(json metadata,
                                                  std::string hty_file_path,
                                                  StrVec_t projected_columns);
void display_one(std::string column_name, Num32Vec &data);
void display_many(StrVec_t column_names, std::vector<Num32Vec> proj_set);
std::optional<UInt32Vec_t> filter_idx(Num32Vec &data, FilterOp op,
                                      Num32 filter_value);
std::optional<Num32Vec> filter_one(Num32Vec &data, FilterOp op,
                                   Num32 filter_value);
std::optional<std::vector<Num32Vec>> filter_many(std::vector<Num32Vec> &data,
                                                 Num32Vec &filtered_col,
                                                 FilterOp op,
                                                 Num32 filter_value);
std::optional<std::vector<Num32Vec>>
project_and_filter(json &metadata, std::string hty_file_path,
                   StrVec_t projected_columns, std::string filtered_column,
                   FilterOp op, Num32 filter_value);
void add_row(json &metadata, std::string hty_file_path,
             std::string modified_hty_file_path,
             std::vector<std::vector<Num32>> rows);

/**
 * @brief Finds the group and column indices for a target column.
 *
 * @param metadata JSON metadata containing group and column information.
 * @param target_col The name of the target column to find.
 * @return A pair of group and column indices if found, or `std::nullopt` if not
 * found.
 */
std::optional<UInt32Pair_t> find_col_idxs(const json &metadata,
                                          std::string target_col)
{
    for (std::uint32_t i = 0; i < metadata["num_groups"]; i++) {
        const json &group = metadata["groups"][i];
        for (std::uint32_t j = 0; j < group["num_columns"]; j++) {
            if (group["columns"][j]["column_name"] == target_col) {
                return (UInt32Pair_t){i, j};
            }
        }
    }
    return std::nullopt;
}

/**
 * @brief Finds the group index and column indices for multiple target columns.
 *
 * @param metadata JSON metadata containing group and column information.
 * @param target_cols A vector of target column names to find.
 * @return A pair of the group index and a vector of column indices, or
 * `std::nullopt` if any column is not found.
 */
std::optional<UInt32ValVecPair_t> find_many_col_idx(const json &metadata,
                                                    StrVec_t target_cols)
{
    std::uint32_t group_idx;
    UInt32Vec_t col_idxs(target_cols.size());
    if (target_cols.size() > 0) {
        std::optional<UInt32Pair_t> idxs_found =
            find_col_idxs(metadata, target_cols[0]);
        if (!idxs_found) {
            std::cerr << "Column '" << target_cols[0] << "' not found\n";
            return std::nullopt;
        }
        group_idx = idxs_found->first;
        col_idxs[0] = idxs_found->second;
    }
    for (std::uint32_t i = 1; i < target_cols.size(); i++) {
        std::optional<std::uint32_t> col_idx = std::nullopt;
        const json &group = metadata["groups"][group_idx];
        for (std::uint32_t j = 0; j < group["num_columns"]; j++) {
            if (group["columns"][j]["column_name"] == target_cols[i]) {
                col_idx = j;
            }
        }
        if (!col_idx) {
            std::cerr << "No Column '" << target_cols[i] << "' in group '"
                      << group_idx << "'\n";
            return std::nullopt;
        }
        col_idxs[i] = *col_idx;
    }
    return (UInt32ValVecPair_t){group_idx, col_idxs};
}

/**
 * @brief Extracts metadata from an HTY file.
 *
 * @param hty_file_path The path to the HTY file.
 * @return A JSON object containing metadata, or `std::nullopt` if extraction
 * fails.
 */
std::optional<json> extract_metadata(std::string hty_file_path)
{
    json header;
    std::ifstream hty(hty_file_path, std::ios::binary);
    if (!hty.is_open()) {
        std::cerr << "Cannot open " << hty_file_path << "\n";
        return std::nullopt;
    }
    try {
        std::uint32_t header_size;
        hty.read((char *)&header_size, sizeof(header_size));
        std::string header_str(header_size, '\0');
        hty.read((char *)&header_str[0], header_size);
        header = json::parse(header_str);
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        hty.close();
        return std::nullopt;
    }
    hty.close();

    auto is_valid_field = [](json &obj, std::string field,
                             bool (json::*type)() const) -> bool {
        if (!obj.contains(field)) {
            std::cerr << "Missing: " << field << "\n";
        } else if (!(obj[field].*type)()) {
            std::cerr << "Incorrect type: " << field << "\n";
        } else {
            return true;
        }
        return false;
    };
    if (!is_valid_field(header, "num_rows", &json::is_number_integer) ||
        !is_valid_field(header, "num_groups", &json::is_number_integer) ||
        !is_valid_field(header, "groups", &json::is_array)) {
        return std::nullopt;
    }
    for (auto &group : header["groups"]) {
        if (!is_valid_field(group, "num_columns", &json::is_number_integer) ||
            !is_valid_field(group, "offset", &json::is_number_integer) ||
            !is_valid_field(group, "columns", &json::is_array)) {
            return std::nullopt;
        }
        for (auto &col : group["columns"]) {
            if (!is_valid_field(col, "column_name", &json::is_string) ||
                !is_valid_field(col, "column_type", &json::is_string)) {
                return std::nullopt;
            }
            try {
                col["column_type"] =
                    num32_type_str_bimap.get_reverse(col["column_type"]);
            } catch (std::out_of_range e) {
                std::cerr << e.what() << "\n";
                return std::nullopt;
            }
        }
    }
    return header;
}

/**
 * @brief Projects a single column from an HTY file.
 *
 * @param metadata JSON metadata containing group and column information.
 * @param hty_file_path The path to the HTY file.
 * @param projected_column The name of the column to project.
 * @return A Num32Vec containing the projected column data, or `std::nullopt` if
 * the operation fails.
 */
std::optional<Num32Vec> project_one(json &metadata, std::string hty_file_path,
                                    std::string projected_column)
{
    std::optional<UInt32Pair_t> idxs =
        find_col_idxs(metadata, projected_column);
    if (!idxs) {
        std::cerr << "Column '" << projected_column << "' not found\n";
        return std::nullopt;
    }

    std::ifstream hty(hty_file_path, std::ios::binary);
    if (!hty.is_open()) {
        std::cerr << "Cannot open " << hty_file_path << "\n";
        return std::nullopt;
    }
    std::vector<Num32> col(metadata["num_rows"]);
    json &group = metadata["groups"][idxs->first];
    Num32Type type = group["columns"][idxs->second]["column_type"];
    std::uint32_t offset =
        group["offset"].get<std::uint32_t>() + idxs->second * sizeof(Num32);
    std::uint32_t row_offset =
        group["num_columns"].get<std::uint32_t>() * sizeof(Num32);

    for (Num32 &cell : col) {
        hty.seekg(offset, std::ios::beg);
        hty.read((char *)&cell, sizeof(cell));
        offset += row_offset;
    }
    hty.close();
    return (Num32Vec){col, type};
}

/**
 * @brief Displays a single column of data.
 *
 * @param column_name The name of the column to display.
 * @param data The Num32Vec containing the column data to display.
 */
void display_one(std::string column_name, Num32Vec &data)
{
    std::cout << column_name << "\n";
    switch (data.type) {
    case Num32Type::I:
        for (Num32 cell : data.vec) {
            std::cout << std::get<std::int32_t>(cell) << "\n";
        }
        break;
    case Num32Type::F:
        for (Num32 cell : data.vec) {
            std::cout << std::get<float>(cell) << "\n";
        }
        break;
    }
}

/**
 * @brief Retrieves a comparison function for a specified filtering operation.
 *
 * @tparam T The type of the values to compare.
 * @param op The filtering operation.
 * @return A function pointer to the comparison function.
 */
template <typename T> static bool (*_get_op(FilterOp op))(T, T)
{
    switch (op) {
    case FilterOp::EQ:
        return [](T a, T b) -> bool { return a == b; };
        break;
    case FilterOp::NEQ:
        return [](T a, T b) -> bool { return a != b; };
        break;
    case FilterOp::L:
        return [](T a, T b) -> bool { return a < b; };
        break;
    case FilterOp::LEQ:
        return [](T a, T b) -> bool { return a <= b; };
        break;
    case FilterOp::G:
        return [](T a, T b) -> bool { return a > b; };
        break;
    case FilterOp::GEQ:
        return [](T a, T b) -> bool { return a >= b; };
        break;
    }
    return [](T a, T b) -> bool { return false; };
}

/**
 * @brief Filters indices of data based on a specified filtering operation and
 * value.
 *
 * @param data The Num32Vec containing the data to filter.
 * @param op The filtering operation.
 * @param filter_value The value to compare against.
 * @return A vector of filtered indices, or `std::nullopt` if filtering fails.
 */
std::optional<UInt32Vec_t> filter_idx(Num32Vec &data, FilterOp op,
                                      Num32 filter_value)
{
    if (data.vec[0].index() != filter_value.index()) {
        std::cerr << "Filter value with type '"
                  << num32_type_str_bimap.get((Num32Type)filter_value.index())
                  << "' does not match type '"
                  << num32_type_str_bimap.get((Num32Type)data.vec[0].index())
                  << "'\n";
        return std::nullopt;
    }
    UInt32Vec_t idxs;
    switch (data.type) {
    case Num32Type::I: {
        auto comp_func = _get_op<std::int32_t>(op);
        for (std::uint32_t i = 0; i < data.vec.size(); i++) {
            if (comp_func(std::get<std::int32_t>(data.vec[i]),
                          std::get<std::int32_t>(filter_value))) {
                idxs.insert(idxs.end(), i);
            }
        }
        break;
    }
    case Num32Type::F: {
        auto comp_func = _get_op<float>(op);
        for (std::uint32_t i = 0; i < data.vec.size(); i++) {
            if (comp_func(std::get<float>(data.vec[i]),
                          std::get<float>(filter_value))) {
                idxs.insert(idxs.end(), i);
            }
        }
        break;
    }
    }
    return idxs;
}

/**
 * @brief Filters data from a Num32Vec based on a specified operation and value.
 *
 * @param data The Num32Vec containing the data to filter.
 * @param op The filtering operation.
 * @param filter_value The value to compare against.
 * @return A Num32Vec containing the filtered data, or `std::nullopt` if
 * filtering fails.
 */
std::optional<Num32Vec> filter_one(Num32Vec &data, FilterOp op,
                                   Num32 filter_value)
{
    std::optional<UInt32Vec_t> idxs = filter_idx(data, op, filter_value);
    if (!idxs) {
        return std::nullopt;
    }
    Num32Vec filtered_data = {.vec = get_vec_subset(data.vec, *idxs),
                              .type = data.type};
    return filtered_data;
}

/**
 * @brief Calculates row offsets for column data in a group.
 *
 * @param col_idxs A vector of column indices.
 * @param total_col The total number of columns in the group.
 * @return A vector of row offsets.
 */
UInt32Vec_t calc_row_offsets(UInt32Vec_t col_idxs, std::uint32_t total_col)
{
    UInt32Vec_t row_offsets(col_idxs.size());
    std::uint32_t i = 0;
    for (; i < col_idxs.size() - 1; i++) {
        row_offsets[i] = (col_idxs[i + 1] - col_idxs[i]) * sizeof(Num32);
    }
    row_offsets[i] = (total_col + col_idxs[0] - col_idxs[i]) * sizeof(Num32);
    return row_offsets;
};

/**
 * @brief Projects multiple columns from an HTY file.
 *
 * @param metadata JSON metadata containing group and column information.
 * @param hty_file_path The path to the HTY file.
 * @param projected_columns A vector of column names to project.
 * @return A vector of Num32Vec objects, each containing the data of a projected
 * column, or `std::nullopt` if the operation fails.
 */
std::optional<std::vector<Num32Vec>> project_many(json metadata,
                                                  std::string hty_file_path,
                                                  StrVec_t projected_columns)
{
    std::optional<UInt32ValVecPair_t> all_idxs =
        find_many_col_idx(metadata, projected_columns);
    if (!all_idxs) {
        return std::nullopt;
    }
    std::uint32_t &group_idx = all_idxs->first;
    UInt32Vec_t &col_idxs = all_idxs->second;
    std::ifstream hty(hty_file_path, std::ios::binary);
    if (!hty.is_open()) {
        std::cerr << "Cannot open " << hty_file_path << "\n";
        return std::nullopt;
    }
    std::vector<Num32Vec> cols(
        col_idxs.size(), {.vec = std::vector<Num32>(metadata["num_rows"])});

    json &group = metadata["groups"][group_idx];
    std::uint32_t total_col = group["num_columns"];
    std::uint32_t start_pos = group["offset"];

    std::uint32_t offset;
    UInt32Vec_t row_offsets;
    offset = start_pos + col_idxs[0] * sizeof(Num32);
    row_offsets = calc_row_offsets(col_idxs, total_col);

    for (std::uint32_t i = 0; i < cols.size(); i++) {
        cols[i].type = group["columns"][col_idxs[i]]["column_type"];
    }
    for (std::uint32_t i = 0; i < cols[0].vec.size(); i++) {
        for (std::uint32_t j = 0; j < cols.size(); j++) {
            hty.seekg(offset, std::ios::beg);
            hty.read((char *)&cols[j].vec[i], sizeof(cols[j].vec[i]));
            offset += row_offsets[j];
        }
    }
    hty.close();
    return cols;
}

/**
 * @brief Displays multiple columns of data.
 *
 * @param column_names A vector of column names to display.
 * @param proj_set A vector of Num32Vec objects containing the column data to
 * display.
 */
void display_many(StrVec_t column_names, std::vector<Num32Vec> proj_set)
{
    if (proj_set.empty()) {
        std::cerr << "projection set empty\n";
        return;
    }
    std::cout << column_names[0];
    for (std::uint32_t i = 1; i < column_names.size(); i++) {
        std::cout << ",\t" << column_names[i];
    }
    std::cout << "\n";
    for (std::uint32_t i = 0; i < proj_set[0].vec.size(); i++) {
        std::cout << get_num32(proj_set[0].vec[i], proj_set[0].type);
        for (std::uint32_t j = 1; j < proj_set.size(); j++) {
            std::cout << ",\t"
                      << get_num32(proj_set[j].vec[i], proj_set[j].type);
        }
        std::cout << "\n";
    }
}

/**
 * @brief Filters multiple columns of data based on a filtered column.
 *
 * @param data The vector of Num32Vec objects containing the columns to filter.
 * @param filtered_col The Num32Vec containing the filtered column.
 * @param op The filtering operation.
 * @param filter_value The value to compare against.
 * @return A vector of filtered Num32Vec objects, or `std::nullopt` if filtering
 * fails.
 */
std::optional<std::vector<Num32Vec>> filter_many(std::vector<Num32Vec> &data,
                                                 Num32Vec &filtered_col,
                                                 FilterOp op,
                                                 Num32 filter_value)
{
    std::optional<UInt32Vec_t> row_idxs =
        filter_idx(filtered_col, op, filter_value);
    if (!row_idxs) {
        return std::nullopt;
    }
    std::vector<Num32Vec> filtered_data(
        data.size(), {.vec = std::vector<Num32>(row_idxs->size())});
    for (std::uint32_t i = 0; i < data.size(); i++) {
        filtered_data[i].type = data[i].type;
        filtered_data[i].vec = get_vec_subset(data[i].vec, *row_idxs);
    }
    return filtered_data;
}

/**
 * @brief Projects and filters data based on specified criteria.
 *
 * @param metadata JSON metadata containing group and column information.
 * @param hty_file_path The path to the HTY file.
 * @param projected_columns A vector of column names to project.
 * @param filtered_column The name of the column to filter.
 * @param op The filtering operation.
 * @param filter_value The value to compare against.
 * @return A vector of filtered and projected Num32Vec objects, or
 * `std::nullopt` if the operation fails.
 */
std::optional<std::vector<Num32Vec>>
project_and_filter(json &metadata, std::string hty_file_path,
                   StrVec_t projected_columns, std::string filtered_column,
                   FilterOp op, Num32 filter_value)
{
    std::uint32_t filtered_idx;
    auto proj_begin = projected_columns.begin();
    auto proj_end = projected_columns.end();
    StrVec_t all_col_names(proj_begin, proj_end);

    /* I group the filter column and project columns so that I can read and
     * filter in one go with project_many and filter_many */
    auto filtered_iter = std::find(proj_begin, proj_end, filtered_column);
    if (filtered_iter != proj_end) {
        filtered_idx = filtered_iter - proj_begin;
    } else {
        all_col_names.insert(all_col_names.end(), filtered_column);
        filtered_idx = projected_columns.size();
    }
    std::optional<std::vector<Num32Vec>> full_set =
        project_many(metadata, hty_file_path, all_col_names);
    if (!full_set) {
        return std::nullopt;
    }
    std::vector<Num32Vec> proj_set =
        slice_vec(*full_set, 0, projected_columns.size());
    std::optional<std::vector<Num32Vec>> filtered_proj =
        filter_many(proj_set, (*full_set)[filtered_idx], op, filter_value);
    return filtered_proj;
}

/**
 * @brief Updates the metadata with a new row count and recalculates offsets.
 *
 * @param metadata The original JSON metadata to update.
 * @param new_row_count The number of new rows to add.
 * @return A JSON object containing the updated metadata.
 */
json update_metadata(json &metadata, std::uint32_t new_row_count)
{
    std::uint32_t group_count = metadata["num_groups"];
    std::uint32_t total_rows =
        metadata["num_rows"].get<std::uint32_t>() + new_row_count;

    auto count_char = [](std::uint32_t num) -> std::uint32_t {
        return std::to_string(num).size();
    };

    /* Update the header except for the offset, since the offset references it
     * own length I decided to calculate it seperately */
    json new_metadata(metadata);
    new_metadata["num_rows"] = total_rows;
    for (json &group : new_metadata["groups"]) {
        group["offset"] = 0;
        for (json &col : group["columns"]) {
            col["column_type"] =
                num32_type_str_bimap.get(col["column_type"].get<Num32Type>());
        }
    }
    std::uint32_t base_size = new_metadata.dump().size() - group_count;

    /* Calculate the base offset without considering their string's lengths */
    std::vector<std::uint32_t> group_offsets(group_count);
    std::uint32_t cum_col_count = 0;
    for (std::uint32_t i = 0; i < group_count; i++) {
        group_offsets[i] = total_rows * cum_col_count * sizeof(Num32) +
                           base_size + sizeof(base_size);
        cum_col_count +=
            new_metadata["groups"][i]["num_columns"].get<std::uint32_t>();
    }

    std::uint32_t total_digit_count = 0;
    for (std::uint32_t offset : group_offsets) {
        total_digit_count += count_char(offset);
    }

    /* Add the length of unaccounted string length, If the string becomes longer
     * repeat the process util no there is no cahnge in length */
    std::list<std::uint32_t> extra_offset(group_count, total_digit_count);
    auto offsets_iter = group_offsets.begin();
    while (!extra_offset.empty()) {
        std::uint32_t new_offset = *offsets_iter + extra_offset.front();
        std::uint32_t diff = count_char(new_offset) - count_char(*offsets_iter);
        extra_offset.pop_front();
        if (diff > 0) {
            extra_offset.insert(extra_offset.end(), group_count - 1, diff);
        }
        *offsets_iter = new_offset;
        offsets_iter++;
        if (offsets_iter == group_offsets.end()) {
            offsets_iter = group_offsets.begin();
        }
    }
    for (std::uint32_t i = 0; i < group_count; i++) {
        new_metadata["groups"][i]["offset"] = group_offsets[i];
    }
    return new_metadata;
}

/**
 * @brief Adds new rows to the HTY file and updates metadata.
 *
 * @param metadata JSON metadata describing the structure and layout of the HTY
 * file.
 * @param hty_file_path Path to the original HTY file.
 * @param modified_hty_file_path Path to the new HTY file with added rows.
 * @param rows A vector of vectors containing the new rows of data. Each inner
 *             vector represents a row, with the column values matching the
 *             metadata-defined types.
 */
void add_row(json &metadata, std::string hty_file_path,
             std::string modified_hty_file_path,
             std::vector<std::vector<Num32>> rows)
{
    /* Check if input rows is valid and seperate it based on the column groups
     */
    std::uint32_t total_cols = 0;
    std::uint32_t group_count = metadata["num_groups"];
    for (json &group : metadata["groups"]) {
        total_cols += group["num_columns"].get<std::uint32_t>();
    }
    if (total_cols != rows[0].size()) {
        std::cerr << "Input must have " << total_cols << " column(s), but has "
                  << rows[0].size() << " column(s)\n";
        return;
    }
    std::vector<std::vector<Num32>> rearraged_rows(group_count);
    for (std::uint32_t i = 0; i < rows.size(); i++) {
        for (std::uint32_t j1 = 0; j1 < group_count; j1++) {
            std::uint32_t col_count = metadata["groups"][j1]["num_columns"];
            for (std::uint32_t j2 = 0; j2 < col_count; j2++) {
                std::uint32_t j = j1 * col_count + j2;

                json &col = metadata["groups"][j1]["columns"][j2];
                Num32Type col_type = col["column_type"];
                if ((Num32Type)rows[i][j].index() != col_type) {
                    std::cerr << "Cell at row " << i << ", column " << j
                              << " must have type "
                              << num32_type_str_bimap.get(col_type)
                              << ", but has type "
                              << num32_type_str_bimap.get(
                                     (Num32Type)rows[i][j].index())
                              << "\n";
                    return;
                }
                rearraged_rows[j1].insert(rearraged_rows[j1].end(), rows[i][j]);
            }
        }
    }
    std::uint32_t old_row_count = metadata["num_rows"];

    /* Read into memory in case the path of the to be modified file is the same
     * as the input */
    std::ifstream hty_in(hty_file_path, std::ios::binary);
    if (!hty_in.is_open()) {
        std::cerr << "Cannot open " << hty_file_path << "\n";
        return;
    }
    std::vector<std::unique_ptr<Num32[]>> old_data(group_count);
    for (std::uint32_t i = 0; i < group_count; i++) {
        json &group = metadata["groups"][i];
        std::uint32_t col_count = group["num_columns"];
        std::uint32_t offset = group["offset"];

        old_data[i] =
            std::unique_ptr<Num32[]>(new Num32[old_row_count * col_count]);
        hty_in.seekg(offset, std::ios::beg);
        hty_in.read((char *)old_data[i].get(),
                    old_row_count * col_count * sizeof(Num32));
    }
    hty_in.close();

    std::ofstream hty_out(modified_hty_file_path, std::ios::binary);
    if (!hty_out.is_open()) {
        std::cerr << "Cannot open " << modified_hty_file_path << "\n";
        return;
    }

    json new_metadata = update_metadata(metadata, rows.size());
    std::string metadata_str = new_metadata.dump();
    std::uint32_t metadata_size = metadata_str.size();

    hty_out.write((char *)&metadata_size, sizeof(metadata_size));
    hty_out.write((char *)&metadata_str[0], metadata_size);
    for (std::uint32_t i = 0; i < group_count; i++) {
        std::uint32_t col_count = metadata["groups"][i]["num_columns"];
        hty_out.write((char *)old_data[i].get(),
                      old_row_count * col_count * sizeof(Num32));
        hty_out.write((char *)&rearraged_rows[i][0],
                      rearraged_rows[i].size() * sizeof(Num32));
    }
    hty_out.close();
}

int main()
{
    /* Sorry I didn't have time get the a program loop working, so I will leave
     * my test code instead */
    std::string hty_file_path = "a.hty", out_path = "b.hty";

    std::optional<json> header = extract_metadata(hty_file_path);
    if (!header) {
        return -1;
    }
    Num32 filter_val;
    FilterOp op;

    std::string col_name = "hi";
    filter_val = 5;
    op = FilterOp::LEQ;
    std::optional<Num32Vec> col_info =
        project_one(*header, hty_file_path, col_name);
    if (col_info) {
        display_one(col_name, *col_info);
        std::cout << "\n";
        std::optional<Num32Vec> filtered_col_info =
            filter_one(*col_info, op, filter_val);
        if (filtered_col_info) {
            display_one(col_name, *filtered_col_info);
            std::cout << "\n";
        }
    }

    StrVec_t col_names = {"1a", "1n", "hi", "1a"};
    filter_val = (float)3.0;
    op = FilterOp::LEQ;
    std::optional<std::vector<Num32Vec>> res_set =
        project_many(*header, hty_file_path, col_names);
    if (res_set) {
        display_many(col_names, *res_set);
        std::cout << "\n";
        std::optional<std::vector<Num32Vec>> filtered_set =
            filter_many(*res_set, (*res_set)[0], op, filter_val);
        if (filtered_set) {
            display_many(col_names, *filtered_set);
            std::cout << "\n";
        }
    }
    col_names = {"1a", "hi", "1n", "hi", "1a"};
    std::string filter_col = "1n";
    filter_val = 6;
    op = FilterOp::EQ;
    std::optional<std::vector<Num32Vec>> pf_set = project_and_filter(
        *header, hty_file_path, col_names, filter_col, op, filter_val);
    if (pf_set) {
        display_many(col_names, *pf_set);
        std::cout << "\n";
    }

    std::vector<Num32Type> col_types;
    for (json &group : (*header)["groups"]) {
        for (json &col : group["columns"]) {
            col_types.insert(col_types.end(),
                             col["column_type"].get<Num32Type>());
        }
    }
#define NEW_ROW_CNT 10
    std::srand(std::time(0));
    std::vector<std::vector<Num32>> new_rows(
        NEW_ROW_CNT, std::vector<Num32>(col_types.size()));
    for (auto &row : new_rows) {
        for (std::uint32_t i = 0; i < col_types.size(); i++) {
            switch (col_types[i]) {
            case Num32Type::I:
                row[i] = (std::int32_t)(std::round(std::rand() % 10));
                std::cout << std::get<std::int32_t>(row[i]) << ", ";
                break;
            case Num32Type::F:
                row[i] = (float)(std::rand() % 10);
                std::cout << "f" << std::get<float>(row[i]) << ", ";
                break;
            }
        }
        std::cout << "\n";
    }
    col_names = {"1n", "hi", "1a", "a"};
    add_row(*header, hty_file_path, out_path, new_rows);
    std::optional<json> new_header = extract_metadata(out_path);
    if (new_header) {
        std::optional<std::vector<Num32Vec>> new_proj =
            project_many(*new_header, out_path, col_names);
        if (new_proj) {
            display_many(col_names, *new_proj);
        }
    }
    return 0;
}

// int main()
// {
//     json metadata;
//     std::string hty_file_path;
//     std::vector<Num32Vec> projected_data;
//     Num32Vec filtered_data;
//     bool metadata_loaded = false;
//     enum Option {
//         LOAD = 1,
//         PROJ_ONE,
//         PROJ_MANY,
//         FILER_ONE,
//         FILTER_MANY,
//         PROJ_AND_FILTER,
//         EXIT
//     };
//     while (true) {
//         std::cout << "\n=== Menu ===\n";
//         std::cout << "1. Load file\n";
//         std::cout << "2. Project one column\n";
//         std::cout << "3. Project multiple columns\n";
//         std::cout << "4. Filter one column\n";
//         std::cout << "5. Filter multiple columns\n";
//         std::cout << "6. Project and filter\n";
//         std::cout << "7. Exit\n";
//         std::cout << "Choose an option: ";

//         int choice;
//         std::cin >> choice;
//         std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//         // Clear input buffer

//         if (choice == EXIT) {
//             std::cout << "Exiting...\n";
//             break;
//         }

//         switch (choice) {
//         case LOAD: {
//             std::cout << "Enter file path: ";
//             std::cin >> hty_file_path;
//             auto result = extract_metadata(hty_file_path);
//             if (result) {
//                 metadata = *result;
//                 metadata_loaded = true;
//                 std::cout << "File loaded successfully.\n";
//             } else {
//                 std::cout << "Failed to load file.\n";
//             }
//             break;
//         }
//         case PROJ_ONE: {
//             if (!metadata_loaded) {
//                 std::cout << "Load a file first.\n";
//                 break;
//             }
//             std::string column_name;
//             std::cout << "Enter column name to project: ";
//             std::cin >> column_name;
//             auto result = project_one(metadata, hty_file_path, column_name);
//             if (result) {
//                 filtered_data = *result;
//                 display_one(column_name, filtered_data);
//             } else {
//                 std::cout << "Projection failed.\n";
//             }
//             break;
//         }
//         case PROJ_MANY: {
//             if (!metadata_loaded) {
//                 std::cout << "Load a file first.\n";
//                 break;
//             }
//             // std::cin.ignore(); // Clear input buffer
//             std::string column_names;
//             std::cout << "Enter column names to project (comma-separated): ";
//             std::getline(std::cin, column_names);

//             std::stringstream ss(column_names);
//             std::string name;
//             StrVec_t columns;
//             while (std::getline(ss, name, ',')) {
//                 columns.push_back(name);
//             }

//             auto result = project_many(metadata, hty_file_path, columns);
//             if (result) {
//                 projected_data = *result;
//                 display_many(columns, projected_data);
//             } else {
//                 std::cout << "Projection failed.\n";
//             }
//             break;
//         }
//         case FILER_ONE: {
//             if (!metadata_loaded) {
//                 std::cout << "Load a file first.\n";
//                 break;
//             }
//             std::string type_str;
//             Num32Type type;
//             Num32 filter_value;
//             int op_choice;
//             FilterOp op;

//             std::cout << "Enter filter type (int, float): ";
//             std::cin >> type_str;
//             try {
//                 type = num32_type_str_bimap.get_reverse(type_str);
//             } catch (std::exception e) {
//                 std::cerr << "Incorrect Type\n";
//                 break;
//             }
//             std::cout << "Enter filter value: ";
//             switch (type) {
//             case Num32Type::I:
//                 std::cin >> std::get<std::int32_t>(filter_value);
//                 break;
//             case Num32Type::F:
//                 std::cin >> std::get<float>(filter_value);
//                 break;
//             }
//             std::cout << "Choose filter operation (0: Equal, 1: Not Equal, 2:
//             "
//                          "Less Than, 3: Less Than or Equal, 4: Greater Than,
//                          " "5: Greater Than or Equal): ";
//             std::cin >> op_choice;
//             try {
//                 op = static_cast<FilterOp>(op_choice);
//             } catch (std::exception e) {
//                 std::cerr << "Incorrect Operation\n";
//                 break;
//             }
//             auto result = filter_one(filtered_data, op, filter_value);
//             if (result) {
//                 filtered_data = *result;
//                 display_one("Filtered Column", filtered_data);
//             } else {
//                 std::cout << "Filter failed.\n";
//             }
//             break;
//         }
//         case FILTER_MANY: {
//             if (!metadata_loaded) {
//                 std::cout << "Load a file first.\n";
//                 break;
//             }
//             if (projected_data.empty()) {
//                 std::cout << "Project data first.\n";
//                 break;
//             }
//             std::string type_str;
//             Num32Type type;
//             Num32 filter_value;
//             int op_choice;
//             FilterOp op;

//             std::cout << "Enter filter type (int, float): ";
//             std::cin >> type_str;
//             try {
//                 type = num32_type_str_bimap.get_reverse(type_str);
//             } catch (std::exception e) {
//                 std::cerr << "Incorrect Type\n";
//                 break;
//             }
//             std::cout << "Enter filter value: ";
//             switch (type) {
//             case Num32Type::I:
//                 std::cin >> std::get<std::int32_t>(filter_value);
//                 break;
//             case Num32Type::F:
//                 std::cin >> std::get<float>(filter_value);
//                 break;
//             }
//             std::cout << "Choose filter operation (0: Equal, 1: Not Equal, 2:
//             "
//                          "Less Than, 3: Less Than or Equal, 4: Greater Than,
//                          " "5: Greater Than or Equal): ";
//             std::cin >> op_choice;
//             try {
//                 op = static_cast<FilterOp>(op_choice);
//             } catch (std::exception e) {
//                 std::cerr << "Incorrect Operation\n";
//                 break;
//             }
//             auto result =
//                 filter_many(projected_data, filtered_data, op, filter_value);
//             if (result) {
//                 projected_data = *result;
//                 display_many({"Filtered Data"}, projected_data);
//             } else {
//                 std::cout << "Filter failed.\n";
//             }
//             break;
//         }
//         case PROJ_AND_FILTER: {
//             if (!metadata_loaded) {
//                 std::cout << "Load a file first.\n";
//                 break;
//             }
//             // std::cin.ignore(); // Clear input buffer
//             std::string column_names;
//             std::cout << "Enter columns to project (comma-separated): ";
//             std::getline(std::cin, column_names);

//             std::stringstream ss(column_names);
//             std::string name;
//             StrVec_t columns;
//             while (std::getline(ss, name, ',')) {
//                 columns.push_back(name);
//             }

//             std::string filter_column;
//             std::string type_str;
//             Num32Type type;
//             Num32 filter_value;
//             int op_choice;
//             FilterOp op;

//             std::cout << "Enter filter column: ";
//             std::cin >> filter_column;
//             std::cout << "Enter filter type (int, float): ";
//             std::cin >> type_str;
//             try {
//                 type = num32_type_str_bimap.get_reverse(type_str);
//             } catch (std::exception e) {
//                 std::cerr << "Incorrect Type\n";
//                 break;
//             }
//             std::cout << "Enter filter value: ";
//             switch (type) {
//             case Num32Type::I:
//                 std::cin >> std::get<std::int32_t>(filter_value);
//                 break;
//             case Num32Type::F:
//                 std::cin >> std::get<float>(filter_value);
//                 break;
//             }
//             std::cout << "Choose filter operation (0: Equal, 1: Not Equal, 2:
//             "
//                          "Less Than, 3: Less Than or Equal, 4: Greater Than,
//                          " "5: Greater Than or Equal): ";
//             std::cin >> op_choice;
//             try {
//                 op = static_cast<FilterOp>(op_choice);
//             } catch (std::exception e) {
//                 std::cerr << "Incorrect Operation\n";
//                 break;
//             }
//             auto result = project_and_filter(metadata, hty_file_path,
//             columns,
//                                              filter_column, op,
//                                              filter_value);
//             if (result) {
//                 projected_data = *result;
//                 display_many(columns, projected_data);
//             } else {
//                 std::cout << "Operation failed.\n";
//             }
//             break;
//         }
//         default:
//             std::cout << "Invalid choice. Please try again.\n";
//         }
//     }
//     return 0;
// }
