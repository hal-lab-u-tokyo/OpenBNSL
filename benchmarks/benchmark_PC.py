from pgmpy.utils import get_example_model
from pgmpy.sampling import BayesianModelSampling

import sys

sys.path.append("/workspace")

from modules.PC_cpp import PCEstimator_cpp
from modules.structural_distance import structural_errors, DAG2CPDAG, PDAG2CPDAG


def load_data(data_type, sample_size):
    model = get_example_model(data_type)
    sampler = BayesianModelSampling(model)
    data = sampler.forward_sample(size=sample_size)
    return model, data


def main():
    model, data = load_data("cancer", 10000)
    print(data)
    best_model = PCEstimator_cpp(data=data, ESS=5.0)
    print(model)
    print(best_model)
    errors = structural_errors(DAG2CPDAG(model), PDAG2CPDAG(best_model))
    print(errors)


if __name__ == "__main__":
    main()
