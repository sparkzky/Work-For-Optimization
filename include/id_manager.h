#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <vector>

struct PersonalInfo {
    // 使用静态char数组
    char name[20];
    char address[100];
    uint32_t age;
    // ... 其他信息
    // 为了模拟一个更大的数据结构，在这里添加一个填充
    char padding[120];
};

class IDManagermentSystem {
private:
    // 根据ID的哈希值（这里简化为取模）获取对应的map
    std::map<uint64_t, uint32_t> &get_map_for_id(uint64_t id) {
        size_t map_index = id % sharded_maps.size();
        return sharded_maps[map_index];
    }
    // 分片的map数组，存储的是具体个人信息的index
    std::vector<std::map<uint64_t, uint32_t>> sharded_maps;
    // 存储所有的个人信息
    std::vector<PersonalInfo> all_person_info;

public:
    static constexpr size_t MAX_ELEMENT_PER_MAP = 9999999;

    IDManagermentSystem(size_t total_elements) {
        size_t num_maps = static_cast<size_t>(std::ceil(
            static_cast<double>(total_elements) / MAX_ELEMENT_PER_MAP));
        if (!num_maps) {
            num_maps = 1;
        }

        sharded_maps.resize(num_maps);

        std::cout << "系统初始化: 总记录数 = " << total_elements << ", 将使用 "
                  << num_maps << " 个map分片." << std::endl;
    }

    void insert_record(uint64_t id, const PersonalInfo &info) {
        all_person_info.push_back(info);
        uint32_t index = static_cast<uint32_t>(all_person_info.size() - 1);

        auto &target_map = get_map_for_id(id);
        target_map[id] = index;
    }

    // 或许这里可以使用std::optional，但为了简化代码，直接返回指针
    const PersonalInfo *find_record(uint64_t id) {
        auto &target_map = get_map_for_id(id);

        auto it = target_map.find(id);
        if (it != target_map.end()) {
            uint32_t index = it->second;
            return &all_person_info[index];
        }

        return nullptr;
    }
};

inline void test_id_manager(const size_t &N, const size_t &M) {
    IDManagermentSystem id_system(N);

    std::vector<uint64_t> ids_to_insert;
    ids_to_insert.reserve(N);

    std::mt19937_64 rng(
        std::chrono::steady_clock::now().time_since_epoch().count());

    uint64_t base_id = 100000000000000000ULL;

    for (size_t i = 0; i < N; ++i) {
        ids_to_insert.push_back(base_id + i);
    }

    std::shuffle(ids_to_insert.begin(), ids_to_insert.end(), rng);

    std::cout << "\n开始批量插入 " << N << " 条记录..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < N; ++i) {
        PersonalInfo info;
        info.age = i % 100;
        snprintf(info.name, sizeof(info.name), "User_%zu", i);
        snprintf(info.address, sizeof(info.address), "Address for user %zu", i);
        id_system.insert_record(ids_to_insert[i], info);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> insert_duration = end_time - start_time;
    double ti = insert_duration.count();
    double bi = static_cast<double>(N) / ti;

    std::cout << "插入完成！" << std::endl;
    std::cout << "总耗时 (Ti): " << ti << " 秒" << std::endl;
    std::cout << "插入带宽 (Bi = N/Ti): " << bi << " 条记录/秒" << std::endl;


    std::vector<uint64_t> ids_to_search;
    ids_to_search.reserve(M);
    for (size_t i = 0; i < M; ++i) {
        ids_to_search.push_back(ids_to_insert[rng() % N]);
    }

    std::cout << "\n开始批量查询 " << M << " 条记录..." << std::endl;
    start_time = std::chrono::high_resolution_clock::now();

    size_t found_count = 0;
    for (size_t i = 0; i < M; ++i) {
        if (id_system.find_record(ids_to_search[i]) != nullptr) {
            ++found_count;
        }
    }

    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> search_duration = end_time - start_time;
    double ts = search_duration.count();
    double bs = static_cast<double>(M) / ts;

    std::cout << "查询完成！" << std::endl;
    std::cout << "成功找到 " << found_count << "/" << M << " 条记录。"
              << std::endl;
    std::cout << "总耗时 (Ts): " << ts << " 秒" << std::endl;
    std::cout << "查询带宽 (Bs = M/Ts): " << bs << " 条记录/秒" << std::endl
              << std::endl;
}
