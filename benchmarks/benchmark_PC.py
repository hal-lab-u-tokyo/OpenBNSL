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
    dataset = ["alarm", "cancer", "earthquake", "survey", "sachs", "child", "alarm"]
    for data_type in dataset:
        print("\n" * 3)
        max_iter = 1
        ave_score = [0 for _ in range(7)]
        for i in range(max_iter):
            model_original, data = load_data(data_type, 10000)
            model_estimated = PCEstimator_cpp(data=data, ESS=5.0)
            errors = structural_errors(
                DAG2CPDAG(model_original), PDAG2CPDAG(model_estimated)
            )
            ave_score = [x + y for x, y in zip(ave_score, errors)]
            print(model_original, model_estimated)
        ave_score = [x / max_iter for x in ave_score]
        print(data_type)
        print(ave_score)
        break


if __name__ == "__main__":
    main()
