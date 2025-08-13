# Documentation

- [./doxygen/](./doxygen/) The backend uses Doxygen; it parses include and outputs to `./doxygen/build/xml`.
- [./sphinx/](./sphinx/) The frontend uses Sphinx; it integrates the backend’s Doxygen-generated XML via the breathe extension.
- Run [../scripts/gen_docs.sh](../scripts/gen_docs.sh) to execute both doxygen and sphinx; the documentation is output to `./sphinx/build/html`.