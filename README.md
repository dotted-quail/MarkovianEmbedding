MarkovianEmbedding
======

## Initial remarks

This is the initial version of the project.

### Running code

Invoke the scripts from a running Julia instance in the package root directory.
Don't forget to activate the project environment by invoking Julia via
```
julia --project=.
```

or by invoking
```
using Pkg
Pkg.activate(".")
```

### Adaptation
- Please check the LICENSE file and replace it by another license if you don't agree with it.
- Please check the `authors` entry in `Project.toml`
- Remove or replace demo scripts, notebooks and `packages/MySubPackage`
- Set up a git repo for sharing the code. All subdirectories and top level files in this initial version should be under version control
- Consider introducing  subdirectories for large simulation results which are not under version control due to their size
- Consider introducing further subdirectories to struture your project
- May be DrWatson.jl can be helpful for you. This project structure is compatible to it with two exceptions:
  - It doesn't rely on `@quickactivate`.  It assumes  that the skills obtained by learning how to work with Julia environments are useful
  - Notebooks are assumed to be Pluto notebooks which can be  version controlled in straigtforward way, unlike Jupyter notebooks

### Initial file structure


```
MyProject/
├── Manifest.toml
├── notebooks
│   └── presentation.jl
├── packages
│   └── MySubPackage
│       ├── Project.toml
│       └── src
│       │    └── MySubPackage.jl
│       └── test
│           └── runtests.jl
├── papers
│   └── publication.tex
├── Project.toml
├── scripts
│   └── example.jl
├── src
│   └── MyProject.jl
└── test
    └── runtests.jl
```
The essential role of these files is as follows:
- `Project.toml`: The project file describes the project on a high level, for example the package/project dependencies and compatibility constraints are listed in the project file. See the [documentation](https://pkgdocs.julialang.org/v1/toml-files/#Project-and-Manifest)
- `Manifest.toml`:  The manifest file is an absolute record of the state of the packages in the environment. It includes exact information about (direct and indirect) dependencies of the project. Given a Project.toml + Manifest.toml pair, it is possible to instantiate the exact same package environment, which is very useful for reproducibility. See the [documentation](https://pkgdocs.julialang.org/v1/toml-files/#Manifest.toml)
- `LICENSE`: License of the project. By default it is the MIT license
- `README.md`: This file
- `src`: Subdirectory for project specific code as part of the MarkovianEmbedding package representing the project.
- `test`: Unit tests for project code in `src`. Could include something from `scripts`, `notebooks`.
- `packages`: Subdirectory containing packages developed locally in this project. At some stage these could be turned into separate Julia packages
- `scripts`, `notebooks`: "Frontend" code for creating project results
- `papers`: "Wer schreibt, der bleibt" - besides of coding, you probably should publish your results...
- `docs`: not in this template, but there belong the sources for the docs created with `Documenter.jl`
