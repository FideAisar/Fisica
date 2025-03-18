import numpy as np
import math

def three_sigma_rule(data):
    mean = np.mean(data).item()
    std = np.std(data).item()
    lower_bound = mean - 3 * std
    upper_bound = mean + 3 * std
    outliers = []
    for i in range(len(data)):
        if data[i] < lower_bound or data[i] > upper_bound:
            outliers.append(data[i])
    return outliers

def bins_builder(data_p3s):
    xmin = min(data_p3s);
    xmax = max(data_p3s)
    n = len(data_p3s)
    bin_width = round((xmax - xmin)/np.sqrt(n))
    xmin_array = [xmin]
    xmax_array = [xmin + bin_width - 1]
    
    while xmax_array[-1] < xmax:
        xmin_array.append(xmax_array[-1] + 1)
        xmax_array.append(xmin_array[-1] + bin_width - 1)
    
    n_bins = len(xmin_array)
    
    xcentr_array = []
    for k in range(n_bins):
        xcentr_array.append((xmax_array[k] + xmin_array[k]) / 2)
    return xmin_array,xmax_array,xcentr_array

def freq_builder(data_p3s,xmin_array,xmax_array,xcentr_array):
    n = len(data_p3s)
    mean = np.mean(data_p3s).item()
    bin_width = round((max(data_p3s) - min(data_p3s))/np.sqrt(n))
    abs_freq_array = [0] * len(xmin_array)
    rel_freq_array = [0] * len(xmin_array)
    freq_density_array = [0] * len(xmin_array)
    xcenr_times_abs_freq_array = [0] * len(xmin_array)
    abs_freq_times_diff_squared_array = [0] * len(xmin_array)



    for i in range(len(xmin_array)):
        for j in range(len(data_p3s)):
            if data_p3s[j] >= xmin_array[i] and data_p3s[j] <= xmax_array[i]:
                abs_freq_array[i] += 1
    for j in range(len(xmin_array)):
        rel_freq_array[j] =  abs_freq_array[j]/n
        freq_density_array[j] =  rel_freq_array[j]/bin_width
        xcenr_times_abs_freq_array[j] = xcentr_array[j]*abs_freq_array[j]
        abs_freq_times_diff_squared_array[j] = abs_freq_array[j]*(((xcentr_array[j]-mean))**2)

    return xcentr_array,abs_freq_array,rel_freq_array,freq_density_array, xcenr_times_abs_freq_array, abs_freq_times_diff_squared_array


def phantom_classes(xcentr_array, abs_freq_array, n, std_dev, mean, bin_width):
    # Calcolo della costante b per la distribuzione normale
    b = (1 / (std_dev * (np.sqrt(2 * np.pi)))).item()
    
    phantom_xcentr_array = xcentr_array.copy()
    phantom_abs_freq_obs = abs_freq_array.copy()
    phantom_arg_exp = [0] * len(xcentr_array)
    phantom_prob_density = [0] * len(xcentr_array)
    phantom_prob = [0] * len(xcentr_array)
    phantom_abs_freq_expe = [0] * len(xcentr_array)
    
    # valori iniziali
    for i in range(len(xcentr_array)):
        phantom_arg_exp[i] = (xcentr_array[i] - mean)**2 / (2 * (std_dev**2))
        phantom_prob_density[i] = b * np.exp(-phantom_arg_exp[i])
        phantom_prob[i] = phantom_prob_density[i] * bin_width
        phantom_abs_freq_expe[i] = n * phantom_prob[i]
    
    goal = n
    shot = np.sum(phantom_abs_freq_expe).item()

    while shot / goal < 0.9995:
        if phantom_prob_density[0] > phantom_prob_density[-1]:
            new_xcentr = phantom_xcentr_array[0] - bin_width
            insert_index = 0
        else:
            new_xcentr = phantom_xcentr_array[-1] + bin_width
            insert_index = len(phantom_xcentr_array)
        
        new_arg_exp = (new_xcentr - mean)**2 / (2 * (std_dev**2))
        new_prob_density = b * np.exp(-new_arg_exp)
        new_prob = new_prob_density * bin_width
        new_abs_freq_expe = new_prob * n
       
        phantom_xcentr_array.insert(insert_index, new_xcentr)
        phantom_abs_freq_obs.insert(insert_index, 0)
        phantom_arg_exp.insert(insert_index, new_arg_exp)
        phantom_prob_density.insert(insert_index, new_prob_density)
        phantom_prob.insert(insert_index, new_prob)
        phantom_abs_freq_expe.insert(insert_index, new_abs_freq_expe)
        
        shot = np.sum(phantom_abs_freq_expe).item()

    return phantom_xcentr_array, phantom_abs_freq_obs, phantom_arg_exp, phantom_prob_density, phantom_prob, phantom_abs_freq_expe

def merger(arrays):
    arrays = [array.copy() for array in arrays]
    radq = [0] * len(arrays[0])
    squared = [0] * len(arrays[0])

    while len(arrays[0]) > 1 and arrays[0][-1]  < 6:
        merged_value = arrays[0][-1] + arrays[0][-2]
        arrays[0][-2] = merged_value
        arrays[0].pop(-1)

        for i in range(1, len(arrays)):
            merged_value = arrays[i][-1] + arrays[i][-2]
            arrays[i][-2] = merged_value
            arrays[i].pop(-1)

    while len(arrays[0]) > 1 and arrays[0][0]  < 6:
        merged_value = arrays[0][0] + arrays[0][1]
        arrays[0][1] = merged_value
        arrays[0].pop(0)

        for i in range(1, len(arrays)):
            merged_value = arrays[i][0] + arrays[i][1]
            arrays[i][1] = merged_value
            arrays[i].pop(0)
    
    for i in range(len(arrays[0])):
        radq[i] = np.sqrt(arrays[1][i]*(1-arrays[2][i]))
        squared[i] = ((arrays[1][i]-arrays[0][i])/radq[i])**2
    
    arrays.append(radq)
    arrays.append(squared)

    return arrays


def read_data(file_path):
    with open(file_path, 'r') as file:
        data = [float(line.strip()) for line in file if line.strip()]
    return data


def main():
    # Path to the data file
    file_path = 'data.txt'

    data = read_data(file_path)
    outliers = three_sigma_rule(data)
    data_p3s = [item for item in data if item not in outliers]

    n = len(data_p3s)
    mean = np.mean(data_p3s).item()
    std_dev = np.std(data_p3s).item()

    
    xmin_array, xmax_array, xcentr_array = bins_builder(data_p3s)
    #print("\nx-min:", xmin_array, "\n")
    #print("x-max:", xmax_array, "\n")
    #print("x-centr:", xcentr_array, "\n")

    bin_width = round((max(data_p3s) - min(data_p3s))/np.sqrt(n))

    xcentr_array,abs_freq_array,rel_freq_array,freq_density_array,xcentr_times_abs_freq_array,abs_freq_times_diff_squared_array = freq_builder(data_p3s,xmin_array,xmax_array,xcentr_array)
    #print("abs freq:", abs_freq_array, "\n")
    #print("rel freq:", rel_freq_array, "\n")
    #print("freq density:", freq_density_array, "\n")
    #print("xcentr * f. ass:",  xcentr_times_abs_freq_array, "\n")
    #print("f. ass * (xcentr - media)^2:", abs_freq_times_diff_squared_array, "\n")

    phantom_xcentr_array,phantom_abs_freq_obs,phantom_arg_exp,phantom_prob_density,phantom_prob,phantom_abs_freq_expe= phantom_classes(xcentr_array,abs_freq_array,n,std_dev,mean,bin_width)
    
    matrix_of_arrays = [phantom_abs_freq_obs,phantom_abs_freq_expe,phantom_prob]
    matrix_of_merged = merger(matrix_of_arrays)
    merged_ok = matrix_of_merged[0]
    merged_ek = matrix_of_merged[1]
    merged_prob = matrix_of_merged[2]
    merged_radq = matrix_of_merged[3]
    merged_squared = matrix_of_merged[4]


    # Prima tabella
    print("\n\n")
    print(f"{'xcentr':^12} {'freq. ass oss':^12} {'arg exp':^12} {'prob_density':^12} {'prob':^12} {'freq. ass. attesa':^12}")
    print("-" * 85)

    for xcentr, abs_freq, rel_freq, freq_dens, molt, squared in zip(
        xcentr_array,abs_freq_array,rel_freq_array,
        freq_density_array,xcentr_times_abs_freq_array,
        abs_freq_times_diff_squared_array    
    ):
        print(f"{xcentr:^12.1f} {abs_freq:^12} {rel_freq:^12.4f} {freq_dens:^12.4f} {molt:^12.4f} {squared:^12.3f}")
    
    # Classi fantasma
    print("\n\n")
    print(f"{'xcentr':^12} {'freq. ass oss':^12} {'arg exp':^12} {'prob_density':^12} {'prob':^12} {'freq. ass. attesa':^12}")
    print("-" * 85)

    for xcentr, abs_freq, arg_exp, prob_dens, prob, abs_freq_expe in zip(
        phantom_xcentr_array, phantom_abs_freq_obs, phantom_arg_exp,
        phantom_prob_density, phantom_prob, phantom_abs_freq_expe
    ):
        print(f"{xcentr:^12.1f} {abs_freq:^12} {arg_exp:^12.4f} {prob_dens:^12.4f} {prob:^12.4f} {abs_freq_expe:^12.3f}")
   
    sum_abs_freq_expe = np.sum(phantom_abs_freq_expe)

    print("-" * 85)
    print(f"{'':^12} {'':^12} {'':^12} {'':^12} {'':^12} {sum_abs_freq_expe:^12.6f}")

    # Tabella merged con chi quadro
    print("\n\n")
    print(f"{'ok':^12} {'ek':^12}  {'prob':^12} {'radq':^12} {'squared':^12}")
    print("-" * 65)
    for ok, ek, prob, radq, squared in zip(
    merged_ok, merged_ek, merged_prob, merged_radq, merged_squared
    ):
        print(f"{ok:^12.0f} {ek:^12.2f} {prob:^12.4f} {radq:^12.4f} {squared:^12.4f}")   
    print("-" * 65)

    chi2 = np.sum(merged_squared)

    print(f"{'':^12} {'':^12}  {'':^12} {'':^12} {chi2:^12.2f}")

    print("\n\n")
    
if __name__ == "__main__":
    main()
