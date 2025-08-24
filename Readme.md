# MatrixLib

A lightweight C++ matrix library implementing dynamic memory management and basic matrix operations, designed as a foundation for more advanced linear algebra functionalities. Currently in development, with core features implemented and many advanced features planned.

---

## Features

In Progress:
- Matrix class with dynamic memory management
- Basic operations: constructors, destructor, and element accessors
- Matrix addition, multiplication, and transpose
- Elementary row operations (swap, scale, add scaled row)
- Property checks (isSquare, isSymmetric, etc.)
- Utility functions (trace, rank)
- Row echelon and reduced row echelon
- Special matrices (identity,zeros,ones,diagonal,scalar,exchange)

Planned features:
- Advanced operations: inverse, LU/QR decomposition, eigenvalues and eigenvectors


---

## About the Project

MatrixLib aims to provide a clean, efficient, and easy-to-use matrix class that supports:

- Basic matrix construction and memory management  
- Common matrix operations such as transpose, addition, and multiplication  
- Elementary row operations (swap, scale, add scaled row)  
- Checks for matrix properties (square, symmetric, identity, etc.)  
- More advanced linear algebra features planned, including inverse, echelon form, LU and QR decompositions, eigenvalues, and eigenvectors

At this stage, only method and function declarations are present. Full implementations, tests, and examples will be added incrementally.

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/DorukSulun/MatrixLib.git
   cd MatrixLib 
   mkdir build && cd build
   cmake ..
   make 
   ```
---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.