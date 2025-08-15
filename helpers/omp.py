import ctypes

class OpenMP:
    """
    Wrapper class for OpenMP functions using ctypes.
    Provides access to core OpenMP threading functionalities.
    """

    def __init__(self, library_name="libgomp.so.1"):
        try:
            self.lib = ctypes.CDLL(library_name)
        except OSError as e:
            raise OSError(f"Failed to load OpenMP library '{library_name}': {e}")

        # Define function prototypes
        self.lib.omp_get_max_threads.restype = ctypes.c_int
        self.lib.omp_get_num_threads.restype = ctypes.c_int
        self.lib.omp_set_num_threads.argtypes = [ctypes.c_int]

    def get_max_threads(self):
        """Returns the maximum number of threads available."""
        return self.lib.omp_get_max_threads()

    def get_num_threads(self):
        """Returns the number of threads in the current parallel region."""
        return self.lib.omp_get_num_threads()

    def set_num_threads(self, num_threads):
        """Sets the number of threads to use in parallel regions."""
        if not isinstance(num_threads, int) or num_threads <= 0:
            raise ValueError("Number of threads must be a positive integer.")
        self.lib.omp_set_num_threads(num_threads)
