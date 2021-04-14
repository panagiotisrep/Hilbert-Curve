# Hilbert Curve

A C++ package for the Hilbert Curve in n dimensions. To use the code, just include the header file in folder 'include'.

To initialize the class:

```
HilbertCurve::HilbertCurve hilbertCurve(n, p);
```

where
- `n` is the number of dimensions
- `p` is the number of iterations


To get the Hilbert value of point (5, 5, 5)
```
std::vector<unsigned long> point{5, 5, 5};
unsigned long hilbertValue = hilbertCurve.hilbertNumberFromPoint(point);
```

To get the point with Hilbert value 100

```
std::vector<unsigned long> point = hilbertCurve.pointFromHilbertNumber(100)
```

You can also create a custom class and functor to sort with respect to Hilbert values (see the example)
```
// create some entities
Entity e1(1, 30, 2, 1000);
Entity e2(2, 32, 6, 1500.5);
Entity e3(3, 40, 15, 780.8);
Entity e4(4, 31, 4, 860.6);
Entity e5(5, 45, 20, 2043.4);

std::vector<Entity *> entities{&e1, &e2, &e3, &e4, &e5};

// sort the entities w.r.t. Hilber value
std::vector<unsigned long> hilbertValues;
hilbertCurve.sortData<Entity, struct myFunctor>(entities, hilbertValues);
```

You can find more documentation in the project's [page](https://panagiotisrep.github.io/Hilbert-Curve).
## Example
There is an example in the folder src. We create a custom class `Entity`, a functor which takes an
instance of that class and maps it to a point in the hypercube lattice in 3D, and then we sort a
vector of `Entity` objects, with respect to their Hilbert value.

To run the example:
```
cmake .
make
./hilbertCurve
```

## Issues

- Let n be the number of dimensions and p the number of iterations. Due to the way we compute the
 binary representation of numbers, it must hold that n*p < size_of(long) * 8. In most computers that would be n*p < 64.
 
