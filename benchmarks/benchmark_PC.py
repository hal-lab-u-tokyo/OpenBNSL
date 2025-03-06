from pgmpy.utils import get_example_model
from pgmpy.sampling import BayesianModelSampling

import sys

sys.path.append("/workspace")

from modules.PC_cpp import PCEstimator_cpp
from modules.RAI_cpp import RAIEstimator_cpp
from modules.structural_distance import structural_errors, PDAG2CPDAG


def load_data(data_type, sample_size):
    model = get_example_model(data_type)
    sampler = BayesianModelSampling(model)
    data = sampler.forward_sample(size=sample_size)
    return model, data


def main():
    dataset = ["andes", "cancer", "earthquake", "survey", "sachs", "child", "alarm"]
    for data_type in dataset:
        print("\n" * 3)
        max_iter = 1
        ave_score = [0 for _ in range(7)]
        for i in range(max_iter):
            model_original, data = load_data(data_type, 100000)
            model_estimated = PCEstimator_cpp(data=data, ESS=5.0)
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
            # erros2 = structural_errors(
            #     PDAG2CPDAG(model_original), PDAG2CPDAG(model_estimated2)
            # )
            print("errors:", errors)
            # print("errors2:", erros2)
            # diff = structural_errors(
            #     PDAG2CPDAG(model_estimated), PDAG2CPDAG(model_estimated2)
            # )
            # print("diff:", diff)
            ave_score = [x + y for x, y in zip(ave_score, errors)]
            print(model_original, model_estimated)
        ave_score = [x / max_iter for x in ave_score]
        print(data_type)
        print(ave_score)
        break


if __name__ == "__main__":
    main()
