#include "../include/id_manager.h"


int main(int argc, char **argv) {
    const size_t N_10e6 = 1000000;
    const size_t M = 1000000;

    test_id_manager(N_10e6, M);
    test_id_manager(N_10e6 * 10, M);
    test_id_manager(N_10e6 * 20, M);
    // test_id_manager(N_10e6 * 100, M);

    return 0;
}
