import numpy as np
from scipy.optimize import curve_fit

class ChlorineDecay:
    def __init__(self, network, emtin_data, method='time', trials=10, accuracy=0.001, delete_extra_points=False):
        self.network = network
        self.emtin_data = emtin_data
        self.method = method
        self.trials = trials
        self.accuracy = accuracy
        self.delete_extra_points = delete_extra_points
        self.current_trial = None
        self.result_log = ''
        self.coeff_dma_result = {}

    def _get_cl_concentration(self):
        # Function to get current chlorine concentration (e.g. at the end of the pipe)
        # Here you can replace this function with your own code to get the chlorine concentration data
        return 1.0  # Just an example, you should replace this with the actual concentration

    def _chlorine_decay(self, t, k):
        c0 = self._get_cl_concentration()
        return c0 * np.exp(-k * t)

    def _run_iterations(self):
        # Step 1: Generate dataframe with initial values
        df_dict = {
            'time': [0.0],
            'conc': [self._get_cl_concentration()]
        }

        # Initialize variables
        interpolated_idx = 0
        demands_list_plt = []

        # 3 ----- run iterations interpolating
        for trial in range(self.trials):
            # Set current trial
            self.current_trial = trial

            # Create dataframe
            df = pd.DataFrame(df_dict)

            # Interpolate NaN values
            interpolated_df = df.interpolate(method=self.method)

            # Fit chlorine decay curve
            popt, pcov = curve_fit(self._chlorine_decay, interpolated_df['time'], interpolated_df['conc'])

            # Calculate current concentration at final time
            t_final = interpolated_df['time'].iloc[-1]
            c_final = self._chlorine_decay(t_final, popt[0])

            # Add results to dataframe
            df_dict['time'].append(t_final + 1)  # Add one more minute
            df_dict['conc'].append(self._get_cl_concentration())

            # Check if concentration has reached target value
            if abs(c_final - self.target_concentration) < self.accuracy:
                # Save results
                self.coeff_dma_result['trial_' + str(trial)] = {}
                self.coeff_dma_result['trial_' + str(trial)]['k'] = popt[0]
                self.coeff_dma_result['trial_' + str(trial)]['final_conc'] = c_final

                # Fill info log
                self.result_log += f'Trial {trial}: K = {popt[0]}, final concentration = {c_final}\n'
                break

            # Delete farthest point from dataframe
            if self.delete_extra_points:
                if interpolated_idx == 1:
                    df_dict['time'].pop(len(df_dict['time']) - 1)
                    df_dict['conc'].pop(len(df_dict['conc']) - 1)
                elif interpolated_idx == 2:
                    df_dict['time'].pop(0)
                    df_dict['conc'].pop(0)
                    interpolated_idx -= 1

            # Now with this dataframe we insert [None, p_f] and interpolate again
            for i in range(len(df_dict['time']) - 1):
                if df_dict['conc'][i] < self.target_concentration < df_dict['conc'][i + 1]:
                    df










