#include <iostream>
#include "src/data.cuh"


int main() {

    Idim3 NodeSize(128, 1, 1);
    Idim3 GridSize(1, 1, 1);
    Idim3 id(0, 0, 0);
    int dim = 3;

    data test(NodeSize, GridSize, id, 0, dim, true);

    std::cout << test.NodeSize.x << ", " << test.NodeSize.y << ", " << test.NodeSize.z << std::endl;
    std::cout << test.map_index(-2, 0, 0, 1) << std::endl;
    std::cout << test.map_index(5, 0, 0, 2) << std::endl;

    test.ToDevice();
}
