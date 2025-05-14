from pgmpy.utils import get_example_model
from pgmpy.sampling import BayesianModelSampling
from pgmpy.models import BayesianNetwork
import pgmpy.metrics
from pgmpy.base import PDAG

import sys

sys.path.append("/workspace")

from modules.PC_cpp import PCEstimator_cpp
from modules.RAI_cpp import RAIEstimator_cpp
from modules.structural_distance import structural_errors, PDAG2CPDAG


def load_data(data_type, sample_size):
    if data_type == "random":
        model = BayesianNetwork.get_random(n_nodes=1000, edge_prob=0.005, n_states=5)
    else:
        model = get_example_model(data_type)
    sampler = BayesianModelSampling(model)
    data = sampler.forward_sample(size=sample_size)
    return model, data

def generate_sample(model, sample_size):
    sampler = BayesianModelSampling(model)
    data = sampler.forward_sample(size=sample_size)
    return data

def main():
    dataset = [
        "munin",
        "cancer",
        "earthquake",
        "survey",
        "sachs",
        "child",
        "alarm",
    ]
    for data_type in dataset:
        print("\n" * 3)
        max_iter = 1
        ave_score = [0 for _ in range(8)]
        errors_all = []
        model_original = get_example_model(data_type)
        for i in range(max_iter):
            print("iter: " + str(i))
            # model_original, data = load_data(data_type, 20000)
            data = generate_sample(model_original, 20000)
            model_estimated, calc_time = PCEstimator_cpp(model=model_original, data=data)
            # model_estimated2, _ = RAIEstimator_cpp(
            #     data=data,
            #     ESS=5.0,
            #     parallel=1,
            #     threshold_DP=0,
            #     search_neighbor=False,
            #     do_orientation_A2=True,
            # )
            errors = structural_errors(
                PDAG2CPDAG(model_original), PDAG2CPDAG(model_estimated)
            )
            errors.append(calc_time)
            # erros2 = structural_errors(
            #     PDAG2CPDAG(model_original), PDAG2CPDAG(model_estimated2)
            # )
            print("errors:", errors)
            errors_all.append(errors)
            # print("errors2:", erros2)
            # diff = structural_errors(
            #     PDAG2CPDAG(model_estimated), PDAG2CPDAG(model_estimated2)
            # )
            # print("diff:", diff)
            ave_score = [x + y for x, y in zip(ave_score, errors)]
            print(model_original, model_estimated)
        ave_score = [x / max_iter for x in ave_score]
        print(data_type)
        print(errors_all)
        print(ave_score)
        break


if __name__ == "__main__":
    main()
