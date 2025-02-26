import numpy as np


def dataframe_wrapper2(data):
    # translate data from pandas dataframe to numpy array and extract column names
    columns = np.array(data.columns.tolist())
    t_data = data.to_numpy()
    # print(columns)

    # get statelist //重い

    statelist = [[] for i in range(len(columns))]
    for i in range(len(columns)):
        for j in range(t_data.shape[0]):
            if t_data[j][i] not in statelist[i]:
                statelist[i].append(t_data[j][i])

    # translate str matrix into int matrix and list of column names into list of int //重い
    data_int = np.zeros(t_data.shape, dtype=np.uint8)
    for i in range(t_data.shape[0]):
        for j in range(len(columns)):
            data_int[i][j] = statelist[j].index(t_data[i][j])
    # print(data_int.dtype)
    # print (t_data.dtype)
    # data_int = openbnsl.str2int_numpy(t_data, statelist)

    n_states = np.zeros(len(columns))
    for i in range(len(columns)):
        n_states[i] = len(statelist[i])
    # print(n_states)
    return columns, data_int, n_states
