/**
 * In this example, we demonstrate how to sort data w.r.t. their Hilbert value. First we create a custom functor,
 * that takes an instance of the class and maps it to the 3D hypercube lattice. Then we call the function sortData
 * to sort the data.
 */
#include "HilbertCurve.h"
#include <iostream>
#include <cmath>
#include <sstream>

namespace hc = HilbertCurve;
typedef hc::HilbertCurve Hilbert;

class Entity {
public:
    Entity(int id, int age, int yearsEmployed, double wage) : id(id), age(age), yearsEmployed(yearsEmployed),
                                                              wage(wage) {};
    int id;
    int age;
    int yearsEmployed;
    double wage;
};

/**
 * Takes an Entity and maps it to the hypercube lattice
 */
struct myFunctor {
    std::vector<unsigned long> operator()(Entity const &entity) {
        return std::vector<unsigned long>{
                static_cast<unsigned long>(entity.age),
                static_cast<unsigned long>(entity.yearsEmployed),
                static_cast<unsigned long>(round(entity.wage))};
    }
} myFunctor;

std::string vecToString(std::vector<unsigned long> vec) {
    std::stringstream ss;
    ss << vec[0] << " " << vec[1] << " " << vec[2];
    return ss.str();
}

int main() {
    Hilbert hilbertCurve(3, 10);

    // create some entities
    Entity e1(1, 30, 2, 1000);
    Entity e2(2, 32, 6, 1500.5);
    Entity e3(3, 40, 15, 780.8);
    Entity e4(4, 31, 4, 860.6);
    Entity e5(5, 45, 20, 2043.4);

    std::vector<Entity *> entities{&e1, &e2, &e3, &e4, &e5};

    // map each entity to a point in the hypercube and compute its Hilbert value
    for (auto e : entities) {
        auto coords = myFunctor(*e);
        std::cout << "Entity " << e->id << " has coordinates " << vecToString(coords)
                  << " and Hilbert value " << hilbertCurve.hilbertNumberFromPoint(coords) << std::endl;
    }

    // sort the entities w.r.t. Hilber value
    std::vector<unsigned long> hilbertValues;
    hilbertCurve.sortData<Entity, struct myFunctor>(entities, myFunctor, hilbertValues);

    int at = 0;
    std::cout << std::endl << "Now Sorted" << std::endl;
    for (auto e : entities) {
        std::cout << "Entity " << e->id << " has Hilbert value " << hilbertValues[at] << std::endl;
        ++at;
    }

    return 0;
}