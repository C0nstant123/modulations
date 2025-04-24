#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include "xil_printf.h"

#define MAX_VALUES 65536
#define NUM_SLOPES 8
#define PI 3.141592653589793

typedef struct {
    double in_phase[MAX_VALUES];
    double quadrature[MAX_VALUES];
    int count;

} ComplexArray;

// Структура для DBSCAN
typedef struct {
    double value;
    int label; // Метка кластера (-1 для шума)
} Point;

ComplexArray read_complex_numbers(const float *realArray, const float *imagArray) {
    ComplexArray complex_array;
    complex_array.count = 0;
    double maxVal=0;
    double absVal=0;

    while (complex_array.count < MAX_VALUES) {
            complex_array.in_phase[complex_array.count] = realArray[complex_array.count];
            complex_array.quadrature[complex_array.count] = imagArray[complex_array.count];

            absVal=sqrt(complex_array.in_phase[complex_array.count] * complex_array.in_phase[complex_array.count] + complex_array.quadrature[complex_array.count] *complex_array.quadrature[complex_array.count]);
            if (absVal > maxVal) {
            	maxVal = absVal;
			}

            complex_array.count++;
    }

    for (size_t i = 0; i < complex_array.count; i++) {
			complex_array.in_phase[i] /= maxVal;
			complex_array.quadrature[i] /= maxVal;
        }


    return complex_array;
}



void calculate_phase(double* i_component, double* q_component, double* t, double* estimated_phase, int len, double fs) {
    for (int i = 0; i < len; i++) {

        // Вычисление фазы с использованием arctan2
        estimated_phase[i] = atan2(q_component[i], i_component[i]);
    }
}
void modulate_baseband(double complex* qpsk_signal, double complex* baseband_signal, int len, double f_c, double* t) {
    for (int i = 1; i < len; i++) {
        double phase = -2 * M_PI * f_c * t[i];
        double complex exp_factor = cos(phase) + I * sin(phase); // exp(-j*2*pi*f_c*t[i])
        baseband_signal[i - 1] = qpsk_signal[i] * exp_factor;
    }
}
double calculate_mean(double* array, int length) {
    double sum = 0.0;

    // Суммирование всех элементов массива
    for (int i = 0; i < length; i++) {
        sum += array[i];
    }

    // Возвращаем среднее значение
    return sum / length;
}


double calculate_std(double* data, int len) {
    double sum = 0.0;
    double mean, variance = 0.0;

    // Вычисление среднего
    for (int i = 0; i < len; i++) {
        sum += data[i];
    }
    mean = sum / len;

    // Вычисление дисперсии
    for (int i = 0; i < len; i++) {
        variance += pow(data[i] - mean, 2);
    }

    // Возвращаем стандартное отклонение
    return sqrt(variance / len);
}

void linear_regression(double* x, double* y, int len, double* slope, double* intercept) {
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;

    for (int i = 0; i < len; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_xx += x[i] * x[i];
    }

    *slope = (len * sum_xy - sum_x * sum_y) / (len * sum_xx - sum_x * sum_x);
    *intercept = (sum_y - (*slope) * sum_x) / len;
}
void sort_array(double* array, int len) {
    for (int i = 0; i < len - 1; i++) {
        for (int j = 0; j < len - i - 1; j++) {
            if (array[j] > array[j + 1]) {
                double temp = array[j];
                array[j] = array[j + 1];
                array[j + 1] = temp;
            }
        }
    }
}
void compute_diff(double* array, double* diff, int len) {
    for (int i = 0; i < len - 1; i++) {
        diff[i] = array[i + 1] - array[i];
    }
}
void compute_diff_int_2(int* array, int* diff, int len) {
    for (int i = 0; i < len - 1; i++) {
        diff[i] = (int)((array[i + 1] - array[i])/2);
    }
}

// //////DBSCAN
double distance(double a, double b) {
    return fabs(a - b); // Поскольку у нас одномерные данные
}

void dbscan(Point* points, int n, double eps, int min_samples) {
    int cluster_id = 0;

    for (int i = 0; i < n; i++) {
        if (points[i].label != 0) { // Если уже обработана
            continue;
        }

        // Найти соседей
        int neighbors[MAX_VALUES];
        int neighbor_count = 0;
        for (int j = 0; j < n; j++) {
            if (distance(points[i].value, points[j].value) <= eps) {
                neighbors[neighbor_count++] = j;
            }
        }

        // Если соседей меньше min_samples, то это шум
        if (neighbor_count < min_samples) {
            points[i].label = -1; // Метка шума
            continue;
        }

        // Создаем новый кластер
        cluster_id++;
        points[i].label = cluster_id;

        // Обрабатываем соседей
        for (int k = 0; k < neighbor_count; k++) {
            int neighbor_idx = neighbors[k];

            if (points[neighbor_idx].label == -1) { // Если был шум
                points[neighbor_idx].label = cluster_id;
            }

            if (points[neighbor_idx].label != 0) { // Если уже обработан
                continue;
            }

            points[neighbor_idx].label = cluster_id;

            // Найти соседей текущего соседа
            int secondary_neighbors[MAX_VALUES];
            int secondary_count = 0;
            for (int m = 0; m < n; m++) {
                if (distance(points[neighbor_idx].value, points[m].value) <= eps) {
                    secondary_neighbors[secondary_count++] = m;
                }
            }

            // Добавляем новых соседей, если их достаточно
            if (secondary_count >= min_samples) {
                for (int m = 0; m < secondary_count; m++) {
                    neighbors[neighbor_count++] = secondary_neighbors[m];
                }
            }
        }
    }
}


int find_max_label(Point* points, int size) {
    int max_label = points[0].label; // Инициализируем первым элементом
    for (int i = 1; i < size; i++) {
        if (points[i].label > max_label) {
            max_label = points[i].label;
        }
    }
    return max_label;
}



void float_to_string(float value, char *str, int precision) {
    int int_part = (int)value;
    float frac_part = value - int_part;
    int frac_int;

    // Преобразуем целую часть
    int n = sprintf(str, "%d", int_part);

    // Преобразуем дробную часть с нужной точностью
    if (precision > 0) {
        str[n] = '.';  // Добавляем точку
        n++;

        // Умножаем дробную часть на 10^precision и округляем
        frac_part *= pow(10, precision);
        frac_int = (int)(frac_part + 0.5);

        // Преобразуем дробную часть в строку
        sprintf(&str[n], "%d", frac_int);
    }
}
int chek_pos_neg(const double *signal, size_t length) {
    int pos = 0;
    int neg = 0;

    // Проверяем первые 100 элементов или длину массива, если она меньше 100
    size_t limit = (length < 100) ? length : 100;

    for (size_t i = 0; i < limit; ++i) {
        if (signal[i] > 0) {
            pos++;
        } else {
            neg++;
        }
    }

    if (pos > 20 && neg > 20) {
        return 2; // И положительных, и отрицательных значений достаточно
    } else if (pos > 20) {
        return 1; // Преобладают положительные значения
    } else if (neg > 20) {
        return 0; // Преобладают отрицательные значения
    }

    return -1; // Никакое из условий не выполнено
}

void process_instantaneous_frequency(
    const double *instantaneous_frequency_pure,
    size_t length,
    double threshold_freq,
    double f_c,
    double *instantaneous_frequency_lora,
    size_t *lora_count) {

    double instantaneous_frequency_neg = 0;
    double instantaneous_frequency_pos = 0;
    int flag_block = 0;
    int flag_normal = 0;

    int direction = (f_c > 0) ? 1 : -1;

    size_t lora_index = 0;

    for (size_t i = 0; i < length; ++i) {
        double freq = instantaneous_frequency_pure[i];

        if (freq > f_c - threshold_freq && freq < f_c + threshold_freq && flag_block == 0) {
            instantaneous_frequency_lora[lora_index++] = freq;
        } else if (freq > f_c - threshold_freq && freq < f_c + threshold_freq && flag_block == 1) {
            flag_normal++;
        }

        if (freq < f_c - threshold_freq && direction > 0 && flag_block == 0) {
            int var_chek = chek_pos_neg(&instantaneous_frequency_pure[i], (i + 100 < length) ? 100 : length - i);
            if (var_chek == 2) {
                flag_block = 1;
//                xil_printf("HAVE NOISE\n");
            } else if (var_chek == 0) {
                flag_block = 1;
//                xil_printf("HAVE ANOTHER SIGNAL\n");
            }
        }

        if (freq > f_c + threshold_freq && direction < 0 && flag_block == 0) {
            int var_chek = chek_pos_neg(&instantaneous_frequency_pure[i], (i + 100 < length) ? 100 : length - i);
            if (var_chek == 2) {
                flag_block = 1;
//                xil_printf("HAVE NOISE\n");
            } else if (var_chek == 1) {
                flag_block = 1;
//                xil_printf("HAVE ANOTHER SIGNAL\n");
            }
        }

        if (instantaneous_frequency_pos > 0 && instantaneous_frequency_neg > 0) {
            flag_block = 1;
        }

        if (flag_normal > 20) {
            flag_block = 0;
            flag_normal = 0;
        }
    }

    *lora_count = lora_index;
}

void moving_average(const double* input, double* output, int size, int window_size) {
    double window_sum = 0;
    int half_window = window_size / 2;
//    xil_printf("window_size %d\n",window_size);
//    xil_printf("size %d\n",size);

    for (int i = 0; i < window_size && i < size; i++) {
        window_sum += input[i];
    }

    for (int i = half_window; i < size; i++) {
        output[i-half_window] = window_sum / window_size;
        if (i - half_window >= 0) {
            window_sum -= input[i - half_window];
        }
        if (i + half_window + 1 < size) {
            window_sum += input[i + half_window + 1];
        }
    }
//    xil_printf("end %d\n",size);
}

void process_signal(double* diff_new_freq_time2, int len, double mean_pila, double std_pila,
                    int* signal_indices_2, int* signal_indices, int* signal_count, double* signal_diff_val, int* diff_val_count) {
    int flag_noise = 0;
    *signal_count = -1; // Счетчик индексов в signal_indices_2
    *diff_val_count = 0; // Счетчик значений в signal_diff_val

    for (int i = 0; i < len; i++) {
        if (diff_new_freq_time2[i] > mean_pila + std_pila || diff_new_freq_time2[i] < mean_pila - std_pila) {
            flag_noise++;
        } else {
            if (flag_noise > 5) {
                (*signal_count)++;
                signal_indices_2[*signal_count] = i - 1;

                if (*signal_count > 0) {
                    int pre_signal_index = signal_indices_2[*signal_count-1]; // Предыдущий индекс
                    int mid_index_offset = (signal_indices_2[*signal_count] - pre_signal_index) / 2;
                    int target_index = pre_signal_index + mid_index_offset;

                    signal_indices[*diff_val_count] = target_index;
                    signal_diff_val[*diff_val_count] = diff_new_freq_time2[target_index];
//                        char str1[20];
//                        sprintf(str1, "%.6f", signal_diff_val[*diff_val_count]);
//                        printf("signal_diff_val111  : %s,\n",str1);
                        printf("signal_indices_ : %d,\n",signal_indices[*diff_val_count]);
                    (*diff_val_count)++;
                }

            }
            flag_noise = 0;
        }

        if (flag_noise == 5) {
            (*signal_count)++;
            signal_indices_2[*signal_count] = i;
//            printf("signal_indices_2 %d  : %d,\n",*signal_count,signal_indices_2[*signal_count]);


        }
    }
}

// Функция для сдвига частоты
void apply_frequency_shift(double complex *signal, double f_shift, int n, double fs) {
    for (int i = 0; i < n; i++) {
        double t = i / fs;
        double complex exponent = cexp(-I * 2 * M_PI * f_shift * t);
        signal[i] *= exponent;
//        signal[i] *= 1000;
    }
}


void butterworth_lowpass_4th_order(double* signal, int LenN,int Fs, int F_high) {
    double wc = tan(PI * F_high / Fs); // Омическая частота
    double k1 = sqrt(2.0) * wc;
    double k2 = wc * wc;

    double a0 = 1 + sqrt(2) * k1 + k2;
    double a1 = 2 * (k2 - 1) / a0;
    double a2 = (1 - sqrt(2) * k1 + k2) / a0;

    double b0 = k2 / a0;
    double b1 = 2 * b0;
    double b2 = b0;

    double filtered_signal[MAX_VALUES];
    filtered_signal[0] = signal[0];
    filtered_signal[1] = signal[1];

    for (int i = 2; i < LenN; i++) {
        filtered_signal[i] = b0 * signal[i] + b1 * signal[i - 1] + b2 * signal[i - 2]
            - a1 * filtered_signal[i - 1] - a2 * filtered_signal[i - 2];

        // Проверка на NaN
        if (isnan(filtered_signal[i])) {
            printf(" NaN in filtred sig %d\n", i);
            return;
        }
    }

    for (int i = 0; i < LenN; i++) {
        signal[i] = filtered_signal[i];
    }
}


void apply_hamming_window(double complex* analytic_signal, int LenN) {
    int k1 = 500;
    int k2 = k1*2;
    for (int i = 0; i < k1; i++) {
        // Наложение окна Хамминга
        double window_value = 0.54 - 0.46 * cos(2 * M_PI * i / (k2 - 1));
        analytic_signal[i] *= window_value;  // Умножаем комплексный сигнал на значение окна
    }
    for (int i = LenN-k1; i < LenN; i++) {
        // Наложение окна Хамминга
        double window_value = 0.54 - 0.46 * cos(2 * M_PI * (i-LenN+2*k1) / (k2 - 1));
        analytic_signal[i] *= window_value;  // Умножаем комплексный сигнал на значение окна
    }
}


void normalize_analytic_signal(double complex *analytic_signal, int length) {
    double mean_real = 0.0, mean_imag = 0.0;
    double std_dev = 0.0;
    double k_norm;

    // Вычисляем среднее значение
    for (int i = 0; i < length; i++) {
        mean_real += creal(analytic_signal[i]);
        mean_imag += cimag(analytic_signal[i]);
    }
    mean_real /= length;
    mean_imag /= length;

    // Центрируем сигнал
    for (int i = 0; i < length; i++) {
        analytic_signal[i] -= mean_real + I * mean_imag;
    }

    // Вычисляем стандартное отклонение
    for (int i = 0; i < length; i++) {
        double real_part = creal(analytic_signal[i]);
        double imag_part = cimag(analytic_signal[i]);
        std_dev += real_part * real_part + imag_part * imag_part;
    }
    std_dev = sqrt(std_dev / length);

    // Нормализация
    k_norm = 1.0 / (std_dev * 2.0);
    for (int i = 0; i < length; i++) {
        analytic_signal[i] *= k_norm;
    }

//    printf("Standard Deviation: %f\n", std_dev);
//    printf("Normalization Factor: %f\n", k_norm);
}



double find_max(double* array, int length) {
    if (length <= 0) return -INFINITY;
    double max_value = array[0];
    for (int i = 1; i < length; i++) {
        if (array[i] > max_value) {
            max_value = array[i];
        }
    }
    return max_value;
}


void process_sf_range(int* sf_range, int sf_count, int bw, int fs, int f_lora, double complex* analytic_signal, int signal_length, double* max_sf ) {

//	  printf("open corel, sf_count=%d,bw=%d,fs=%d,f_llora=%d, signal_length=%d; max_sf=%f    \n",sf_count,bw,fs,f_lora,signal_length,max_sf[0]);
//    for (int i=0;i<6;i++){
//        printf("analytic_signal i = %d",sf_range[i]);
//    }
//    printf("\n");

    for (int sf_index = 0; sf_index < sf_count; sf_index++) {
        int sf = sf_range[sf_index];
        double Ts = pow(2, sf) / bw;  // Период символа
        double k = bw / Ts;  // Крутизна чирпа
;
        // Временная ось для одного символа
        int Ns = (int)(Ts * fs);  // Количество отсчетов для одного символа
        double* t_sym = (double*)malloc(Ns * sizeof(double));
        double complex* lora_chirp = (double complex*)malloc(Ns * sizeof(double complex));


        double  lora_chirp3 [MAX_VALUES];

        for (int i = 0; i < Ns; i++) {
            t_sym[i] = (double)i / fs;
            lora_chirp[i] = cexp(I * 2 * M_PI * (f_lora * t_sym[i] + 0.5 * k * t_sym[i] * t_sym[i]));

        }

//        printf("\n", Ns);

        double complex* lora_chirp2 = (double complex*)malloc(signal_length * sizeof(double complex));
        for (int i = 0; i < signal_length; i++) {
            lora_chirp2[i] = lora_chirp[i % Ns];
        }

        int lag_step = 1000;
        int lag_count = (2 * signal_length - 1) / lag_step + 1;
        double* correlation = (double*)malloc(lag_count * sizeof(double));

        for (int i = 0, lag = -signal_length + 1; i < lag_count; i++, lag += lag_step) {
            double complex corr_value = 0.0;
            if (lag < 0) {
                for (int j = 0; j < signal_length + lag; j++) {
                    corr_value += conj(analytic_signal[j]) * lora_chirp2[j - lag];
                }
            } else if (lag > 0) {
                for (int j = 0; j < signal_length - lag; j++) {
                    corr_value += conj(analytic_signal[j + lag]) * lora_chirp2[j];
                }
            } else {
                for (int j = 0; j < signal_length; j++) {
                    corr_value += conj(analytic_signal[j]) * lora_chirp2[j];
                }
            }
            correlation[i] = creal(corr_value);

        }
//        for (int i=0;i<20;i++)
//            printf("correlation i = %f, ",correlation[i]);
//        printf("\n");

        max_sf[sf - 7] = find_max(correlation, lag_count);


//        for (int i=0;i<6;i++)
//            printf("max_sf i = %f",max_sf[i]);
//        printf("\n");
        free(t_sym);
        free(lora_chirp);
        free(lora_chirp2);
        free(correlation);
    }
}




int CHEK_LORA_MOD_2(const float *realArray, const float *imagArray,double f_c_in, double fs_in)
{
    ComplexArray result = read_complex_numbers(realArray, imagArray);

	xil_printf("Open CHEK_LORA_MOD  2\r\n");

	double t[MAX_VALUES];

	double fs = fs_in;
	double f_c = f_c_in;
	int len_array = MAX_VALUES;

    double complex analytic_signal[MAX_VALUES];
    double real_signal[MAX_VALUES];
	double imag_signal[MAX_VALUES];
	int N = MAX_VALUES;



	for (int i = 0; i < result.count; i++) {
		analytic_signal[i] = result.in_phase[i] + I * result.quadrature[i];  // Используем I для мнимой части
	}

//	for (int i=0;i<26;i++){
//	        printf("analytic_signal i = %f\n",creal(analytic_signal[i]));
//	    }
//	for (int i=0;i<26;i++){
//		        printf("result.in_phase i = %f\n",result.in_phase[i]);
//		    }

	apply_frequency_shift(analytic_signal, f_c, N, fs);

	for (int i = 0; i < result.count; i++) {
		real_signal[i] = creal(analytic_signal[i]);
		imag_signal[i] = cimag(analytic_signal[i]);
	}
	int  F_high = 2000000;
	butterworth_lowpass_4th_order(real_signal, N,fs,F_high);
	butterworth_lowpass_4th_order(real_signal, N,fs,F_high);

	butterworth_lowpass_4th_order(imag_signal, N,fs,F_high);
	butterworth_lowpass_4th_order(imag_signal, N,fs,F_high);

	for (int i = 0; i < result.count; i++) {
		analytic_signal[i] = real_signal[i] + I * imag_signal[i];  // Используем I для мнимой части
	}

	apply_frequency_shift(analytic_signal, -f_c, N, fs);

	apply_hamming_window(analytic_signal, N) ;
	normalize_analytic_signal(analytic_signal, N);

	int sf_range[6] = {7,8,9,10,11,12};
	const int sf_count = 6;
	int bw = 512e3;
	double f_lora[MAX_VALUES];
	double max_sf[sf_count];

	process_sf_range(sf_range,  sf_count, bw,  fs,  f_c, analytic_signal, N, max_sf);


	double max_sf_high;
	max_sf_high = find_max(max_sf,6);

//	printf("\n\n max_sf = %f \n\n",max_sf_high);
	if (max_sf_high>2200){
		xil_printf("There is a LORA modulation!  2\r\n");
	}
	else{
		xil_printf("NOT LORA!  \r\n");
	}





}



int CHEK_LORA_MOD(const float *realArray, const float *imagArray,double f_c_in, double fs_in)
{
//	for (int i = 0; i < 30; i++) {
//		char str[20];
//		float_to_string(realArray[i], str, 6);
//		xil_printf("realArray[%d]: %s\n",i, str);
//
//	}

    ComplexArray result = read_complex_numbers(realArray, imagArray);

	xil_printf("Open CHEK_LORA_MOD  \r\n");

    int length = sizeof(result.in_phase) / sizeof(result.in_phase[0]);
    //

    //

	//
	double t[MAX_VALUES];

	double fs = fs_in;
	double f_c = f_c_in;

	int len_array = MAX_VALUES;

	double complex qpsk_signal[MAX_VALUES];
	double complex baseband_signal[MAX_VALUES - 1];
	double phase_base[MAX_VALUES];
	double phase_base2[MAX_VALUES];
	double phase_base_diff[MAX_VALUES - 1];
	double instantaneous_frequency[MAX_VALUES - 1];  // РњР°СЃСЃРёРІ СЂР°Р·РЅРѕСЃС‚РµР№ С„Р°Р·
    double phase_base_diff_std;
    int times_change_phase[MAX_VALUES];  // РњР°СЃСЃРёРІ РёРЅРґРµРєСЃРѕРІ РёР·РјРµРЅРµРЅРёР№ С„Р°Р·С‹
    int times_change_phase_diff[MAX_VALUES];  // РњР°СЃСЃРёРІ РёРЅРґРµРєСЃРѕРІ РёР·РјРµРЅРµРЅРёР№ С„Р°Р·С‹
    double last_times_phase = 0;
    int times_change_phase_count = 0;  // РЎС‡С‘С‚С‡РёРє РєРѕР»РёС‡РµСЃС‚РІР° РёР·РјРµРЅРµРЅРёР№ С„Р°Р·С‹
    double instantaneous_frequency_lora[MAX_VALUES]; // Р РµР·СѓР»СЊС‚Р°С‚
    size_t lora_count = 0;
    double estimated_phase[MAX_VALUES];


    double count_phase[MAX_VALUES];
    double count_phase_diff[MAX_VALUES];


	// Р›РћР Рђ ..............................................................
	float freqs[MAX_VALUES];
	float spectrum_magnitude[MAX_VALUES];

	for (int i = 0; i < len_array; i++) {
		 // Р’С‹С‡РёСЃР»РµРЅРёРµ РІСЂРµРјРµРЅРё РґР»СЏ РєР°Р¶РґРѕРіРѕ РёРЅРґРµРєСЃР°
		t[i] = i / fs;

		// Р’С‹С‡РёСЃР»РµРЅРёРµ С„Р°Р·С‹ СЃ РёСЃРїРѕР»СЊР·РѕРІР°РЅРёРµРј arctan2
		estimated_phase[i] = atan2(result.quadrature[i], result.in_phase[i]);
	}
//    for (int i = 0; i < result.count && i < 10; i++) {
//        xil_printf("t: %.12f,\n", t[i]);
//    }

	// Р’С‹С‡РёСЃР»РµРЅРёРµ РєРѕРјРїР»РµРєСЃРЅРѕРіРѕ СЃРёРіРЅР°Р»Р°
	for (int i = 0; i < result.count ; i++) {
		qpsk_signal[i] = result.in_phase[i] + I * result.quadrature[i];
	}
	// Р—РґРµСЃСЊ РЅСѓР¶РЅРѕ Р·Р°РїРѕР»РЅРёС‚СЊ analytic_signal


	modulate_baseband(qpsk_signal, baseband_signal, result.count, f_c, t);



	for (int i = 0; i < len_array; i++) {
		// Р’С‹С‡РёСЃР»РµРЅРёРµ С„Р°Р·С‹ РІ СЂР°РґРёР°РЅР°С… Рё РїСЂРµРѕР±СЂР°Р·РѕРІР°РЅРёРµ РІ РіСЂР°РґСѓСЃС‹
		phase_base[i] = (atan2(cimag(baseband_signal[i]), creal(baseband_signal[i])) * 180.0) / M_PI;
		phase_base2[i] = (atan2(result.quadrature[i], result.in_phase[i])) ;
	}

	double add_phase=0;
	for (size_t i = 1; i < len_array; ++i) {
		phase_base2[i] += 2 * M_PI * add_phase;
		double delta = phase_base2[i] - phase_base2[i - 1];
		if (delta < 0) {
			add_phase += 1;
			phase_base2[i] += 2 * M_PI;
		}
	}


	xil_printf("len_array 1 = %d. \n",len_array);
	char str1[20];
	float_to_string(fs, str1, 6);
	xil_printf("fs = %s. ",str1);

	for (int i = 1; i < len_array; i++) {
		phase_base_diff[i-1] = phase_base[i] - phase_base[i-1];
		instantaneous_frequency[i - 1] = (phase_base2[i] - phase_base2[i - 1]) * fs / (2.0 * M_PI);
	}
	for (int i; i<20;i++){
		char str1[20];
		float_to_string(instantaneous_frequency[i], str1, 6);
		xil_printf("instantaneous_frequency = %s. ",str1);
	}


	int window_size_1 = 500;
	double instantaneous_frequency_moving[len_array];  //
	double moving_avg[len_array];  //
	moving_average(instantaneous_frequency, instantaneous_frequency_moving, len_array, window_size_1);

	len_array=len_array-window_size_1;


	process_instantaneous_frequency(
			instantaneous_frequency_moving,
			len_array,
			1e6,
			f_c,
			instantaneous_frequency_lora,
			&lora_count
		);



	xil_printf("\n");
	for (int i; i<20;i++){
		char str2[20];
		float_to_string(instantaneous_frequency_moving[i], str2, 6);
		xil_printf("instantaneous_frequency_moving = %s.  ",str2);
	}
	xil_printf("\n");

//	xil_printf("instantaneous_frequency_lora len: %zu\n", lora_count);
	double moving_std[lora_count];  // РњР°СЃСЃРёРІ РёРЅРґРµРєСЃРѕРІ РёР·РјРµРЅРµРЅРёР№ С„Р°Р·С‹

	double moving_average_diff[lora_count-1];  //


	for (int p=0;p<2;p++){
	if (p==0){
		compute_diff(instantaneous_frequency_lora, moving_average_diff, lora_count-1);
	}
	else{
		int window_size_3 = 30;
		moving_average(moving_avg, moving_avg, lora_count, window_size_3);
		compute_diff(moving_avg, moving_average_diff, lora_count-window_size_1-window_size_3-1);
		int window_size_1 = window_size_1+window_size_3;
	}


		int window_size_2 = 50;
		double diff_new_freq_time2[lora_count-1];  //
		moving_average(moving_average_diff, diff_new_freq_time2, lora_count-window_size_1-1, window_size_2);
		double mean_pila;
		double std_pila;
		int signal_count;
		int diff_val_count;
		int signal_indices_2[ lora_count-window_size_1-window_size_2-1];
		int signal_indices[ lora_count-window_size_1-window_size_2-1];
		double signal_diff_val[ lora_count-window_size_1-window_size_2-1];
		mean_pila = calculate_mean(diff_new_freq_time2,lora_count-window_size_1-window_size_2-1);
		std_pila = calculate_std(diff_new_freq_time2,lora_count-window_size_1-window_size_2-1);

		process_signal(diff_new_freq_time2, lora_count-window_size_1-window_size_2-1, mean_pila, std_pila, signal_indices_2,signal_indices, &signal_count, signal_diff_val, &diff_val_count);

		int count_signal_indices=diff_val_count;
		int signal_indices_2_1[count_signal_indices];
		int signal_indices_2_2[count_signal_indices];
		int signal_indices_3[count_signal_indices];
		int count_signal_indices_2_1 = 0, count_signal_indices_2_2 = 0;

		xil_printf("count_signal_indices = %d. lora_count = %d. window_size_2 = %d \n",count_signal_indices,lora_count,window_size_2);
		if (count_signal_indices>2 && count_signal_indices<5000){
			int min_samples = 2;
			int max1 = 0;
			int max2 = 0;
			for (int i = 0; i < count_signal_indices; i++) {
				signal_indices_3[i] = signal_indices[i];
				if (i % 2 == 0) {
					signal_indices_2_1[count_signal_indices_2_1++] = signal_indices[i];
				} else {
					signal_indices_2_2[count_signal_indices_2_2++] = signal_indices[i];
				}
			}
			if (count_signal_indices>5){
				int signal_indices_2_1_diff[count_signal_indices];
				int signal_indices_2_2_diff[count_signal_indices];
				compute_diff(signal_indices_2_1, signal_indices_2_1_diff, count_signal_indices_2_1);
				compute_diff(signal_indices_2_2, signal_indices_2_2_diff, count_signal_indices_2_2);

				sort_array(signal_indices_2_1_diff, count_signal_indices_2_1-1);
				sort_array(signal_indices_2_2_diff, count_signal_indices_2_2-1);

				double eps1 = signal_indices_2_1_diff[count_signal_indices_2_1 - 2] * 0.10; // max = РїРѕСЃР»РµРґРЅРёР№ СЌР»РµРјРµРЅС‚ РїРѕСЃР»Рµ СЃРѕСЂС‚РёСЂРѕРІРєРё
				double eps2 = signal_indices_2_2_diff[count_signal_indices_2_2 - 2] * 0.10;

				char str1[20];
				float_to_string(eps1, str1, 6);
				char str2[20];
				float_to_string(eps2, str2, 6);

				xil_printf("eps1 = %s. eps2 = %s.  \n",str1,str2);
				for (int i; i<count_signal_indices;i++){
					char str4[20];
					float_to_string(signal_indices_3[i], str4, 6);
					xil_printf("signal_indices_3 = %s. \n",str4);
				}


				int labels_ind1 [count_signal_indices_2_1-1];
				int labels_ind2 [count_signal_indices_2_2-1];

				Point points_int1[MAX_VALUES];
				for (int i = 0; i < count_signal_indices_2_1-1; i++) {
					points_int1[i].value = signal_indices_2_1_diff[i];
					points_int1[i].label = 0; // Р’СЃРµ С‚РѕС‡РєРё РёР·РЅР°С‡Р°Р»СЊРЅРѕ РЅРµ РєР»Р°СЃСЃРёС„РёС†РёСЂРѕРІР°РЅС‹
				}
				Point points_int2[MAX_VALUES];
				for (int i = 0; i < count_signal_indices_2_2-1; i++) {
					points_int2[i].value = signal_indices_2_2_diff[i];
					points_int2[i].label = 0; // Р’СЃРµ С‚РѕС‡РєРё РёР·РЅР°С‡Р°Р»СЊРЅРѕ РЅРµ РєР»Р°СЃСЃРёС„РёС†РёСЂРѕРІР°РЅС‹
				}
				dbscan(points_int1, count_signal_indices_2_1 - 1, eps1, min_samples);
				dbscan(points_int2, count_signal_indices_2_2 - 1, eps2, min_samples);


				for (int i = 1; i < signal_indices_2_1 - 1; ++i) {
				   if (points_int1[max1].label < points_int1[i].label) {
					   max1 = i;
				   }
				}
				for (int i = 1; i < signal_indices_2_2 - 1; ++i) {
				   if (points_int2[max2].label < points_int2[i].label) {
					   max2 = i;
				   }
				}
				if ((-1 < points_int1[max1].label <= 2) || (-1 < points_int2[max2].label <= 2)){
					xil_printf("Lora modulations\n");
					break;
				}
			}

			int signal_indices_diff[count_signal_indices];
			compute_diff(signal_indices_3, signal_indices_diff, count_signal_indices);
			sort_array(signal_indices_diff, count_signal_indices-1);
			double eps3 = signal_indices_diff[count_signal_indices - 2] * 0.10;
			int labels_ind3 [count_signal_indices-1];
			Point points_int3[MAX_VALUES];
			for (int i = 0; i < count_signal_indices-1; i++) {
				points_int3[i].value = signal_indices_diff[i];
				points_int3[i].label = 0; // Р’СЃРµ С‚РѕС‡РєРё РёР·РЅР°С‡Р°Р»СЊРЅРѕ РЅРµ РєР»Р°СЃСЃРёС„РёС†РёСЂРѕРІР°РЅС‹
			}
			dbscan(points_int3, count_signal_indices - 1, eps3, min_samples);
			int max3 = 0;
			for (int i = 1; i < count_signal_indices - 1; ++i) {
				if (points_int3[max3].label < points_int3[i].label) {
				   max3 = i;
				}
			}
			xil_printf("max1 = %d.max2 = %d.max3 = %d \n",max1,max2,max3);

			if (-1 < points_int3[max3].label <= 2){
				xil_printf("Lora modulations\n");
				break;
			}


		}
		else{
			xil_printf("Not Lora\n");
		}
	}
}

int CHEK_FASE_MOD(const float *realArray, const float *imagArray,double f_c_in, double fs_in)
{

	ComplexArray result = read_complex_numbers(realArray, imagArray);

	xil_printf("Open CHEK_FASE_MOD  \r\n");

    // Пример вывода первых 20 значений
//    for (int i = 0; i < result.count && i < 10; i++) {
//    	char str1[20];
//    	char str2[20];
//		sprintf(str1, "%f", result.in_phase[i]);
//		sprintf(str2, "%f", result.quadrature[i]);
//
////        printf("In-phase: %s, Quadrature: %f\n", str1, str2);
//    }
    int length = sizeof(result.in_phase) / sizeof(result.in_phase[0]);
    // Вычисление времени для каждого индекса

    // Обьявление переменных
    double t[MAX_VALUES];

    double fs = fs_in;//частота дескретизации
    double f_c = f_c_in;//Центральная частота сигнала

    int len_array = MAX_VALUES;

    double complex qpsk_signal[MAX_VALUES];
    double complex baseband_signal[MAX_VALUES - 1];
    double phase_base[MAX_VALUES];
    double phase_base_diff[MAX_VALUES - 1];  // Массив разностей фаз
    double phase_base_diff_std;
    int times_change_phase[MAX_VALUES];  // Массив индексов изменений фазы
    int times_change_phase_diff[MAX_VALUES];  // Массив индексов изменений фазы
    double last_times_phase = 0;
    int times_change_phase_count = 0;  // Счётчик количества изменений фазы

    double estimated_phase[MAX_VALUES];
    // Выделяем данные для линейной регрессии
    int slope1_len = 0;
    double intercept;
    double x[MAX_VALUES], y[MAX_VALUES],slope[NUM_SLOPES], slope_2[NUM_SLOPES];;
    int start_slope1_idx = 0;
    int end_slope1_idx = 0;
    int len_slope1_inx = 0;
    int t_lin[500];
    double slope_list_diff[NUM_SLOPES - 1];     // Массив разностей
    double slope_list_diff_2[NUM_SLOPES - 1];   // Второй массив разностей (копия)
    double std_slope_list_diff;
    int shift_freq = 0;
    double slope_real_2;

    double count_phase[MAX_VALUES];
    double count_phase_diff[MAX_VALUES];


    for (int i = 0; i < len_array; i++) {
         // Вычисление времени для каждого индекса
        t[i] = i / fs;

        // Вычисление фазы с использованием arctan2
        estimated_phase[i] = atan2(result.quadrature[i], result.in_phase[i]);
    }

    // Вычисление комплексного сигнала
    for (int i = 0; i < result.count ; i++) {
        qpsk_signal[i] = result.in_phase[i] + I * result.quadrature[i];
    }

    modulate_baseband(qpsk_signal, baseband_signal, result.count, f_c, t);


    for (int i = 0; i < len_array; i++) {
        // Вычисление фазы в радианах и преобразование в градусы
        phase_base[i] = (atan2(cimag(baseband_signal[i]), creal(baseband_signal[i])) * 180.0) / M_PI;
    }
    for (int i = 0; i < len_array - 1; i++) {
        phase_base_diff[i] = phase_base[i + 1] - phase_base[i];
    }


    phase_base_diff_std = calculate_std(phase_base_diff, len_array - 1) / 2;
//    printf("phase diff std= %f",phase_base_diff_std);

    for (int i = 0; i < len_array - 1; i++) {
        if ((phase_base_diff[i] > phase_base_diff_std || phase_base_diff[i] < -phase_base_diff_std) && (t[i] - (double)last_times_phase) > 4 * t[1]) {
            times_change_phase[times_change_phase_count++] = i;
            last_times_phase = t[i];
        }

        if ((phase_base_diff[i] > phase_base_diff_std || phase_base_diff[i] < -phase_base_diff_std) && (t[i] - last_times_phase) < 4 * t[1]) {
            last_times_phase = t[i];
        }
    }
    int times_change_len=times_change_phase_count;
    compute_diff_int_2(times_change_phase, times_change_phase_diff, times_change_len);
	if (times_change_len > 30 ){
//		for (int i = 0; i < result.count && i < 10; i++) {
//			xil_printf("compute_diff_int_2 %d  \r\n",times_change_phase_diff[i]);
//			xil_printf("times_change_len %d  \r\n",times_change_len);
//		}
		int count_break = 0;
		int count_break2 = 0;
		// /////////////////////////////////По углу наклона фазы вычисляем неточность фазы

		while (1){
			count_break++;
			count_break2++;

	//        if ((count_break2 % 3) ==1) {
	//			char str1[20];
	//			sprintf(str1, "%f", slope_real_2);  // 2 знака после запятой
	//			xil_printf("slope_real_2]: %s\n", str1);
	//			xil_printf("shift_freq]: %d\n", shift_freq);
	//
	//		}

			if (count_break>=15){
				xil_printf("NOT PSK MODULATIONS!\n");
				break;
			}


			modulate_baseband(qpsk_signal, baseband_signal, result.count, f_c - shift_freq, t);

			for (int i = 0; i < len_array; i++) {
				// Вычисление фазы в радианах и преобразование в градусы
				phase_base[i] = (atan2(cimag(baseband_signal[i]), creal(baseband_signal[i])) * 180.0) / M_PI;
			}

			for (int j = 0; j<NUM_SLOPES; j++){
				slope1_len = (int)((times_change_phase[j+1] - times_change_phase[j]) * 0.2);
				start_slope1_idx = times_change_phase[j] + slope1_len;
				end_slope1_idx = times_change_phase[j+1] - slope1_len;
				len_slope1_inx = end_slope1_idx - start_slope1_idx;
				for (int i = 0; i < 500; i++) {
						t_lin[i] = i + 1;
					}
				int k=0;
				for (int i = start_slope1_idx; i <= end_slope1_idx; i++) {
						y[k] = phase_base[i];
						x[k] = k+1;
						k=k+1;
					}
				linear_regression(x, y, len_slope1_inx, &slope[j], &intercept);
	//            printf("Slope: %f, Intercept: %f\n", slope[j], intercept);
			}
			sort_array(slope, NUM_SLOPES);
			compute_diff(slope, slope_list_diff, NUM_SLOPES);
			for (int i = 0; i < NUM_SLOPES - 1; i++) {
				slope_list_diff_2[i] = slope_list_diff[i];
	//            printf("\n slope_list_diff_2: %f\n", slope_list_diff_2[i]);
			}
			sort_array(slope_list_diff_2, NUM_SLOPES - 1);
			std_slope_list_diff = slope_list_diff_2[1] * 10;
	//        printf("\n std_slope_list_diff: %f\n", std_slope_list_diff);

			int slope_list_2_len = 0;
			for (int i = 0; i < NUM_SLOPES - 1; i++) {
			   if (slope_list_diff[i] < std_slope_list_diff) {
				   slope_2[slope_list_2_len++] = slope[i];
			   }
			   if (slope_list_diff[i] > std_slope_list_diff && slope_list_2_len != 0) {
				   slope_2[slope_list_2_len++] = slope[i];
				   break;
			   }
			   if (i == NUM_SLOPES - 2) {
				   slope_2[slope_list_2_len++] = slope[i + 1];
			   }
			}
	//        printf("Filtered slope_list_2 values:\n");
	//        for (int i = 0; i < slope_list_2_len; i++) {
	//            printf("%f ", slope_2[i]);
	//        }
			slope_list_2_len = sizeof(slope_2) / sizeof(slope_2[0]);
			slope_real_2 = calculate_mean(slope_2, slope_list_2_len);

			if (slope_real_2 < 0) {
				if (slope_real_2 > -0.01) {
					shift_freq += 500;
					if (slope_real_2 > -0.005) {
						break;  // Выход из цикла
					}
				} else {
					shift_freq += (int)(-slope_real_2/0.02*1000);
				}
			}
			if (slope_real_2 > 0) {
				if (slope_real_2 < 0.01) {
					shift_freq -=500;
					if (slope_real_2 < 0.005) {
						break;  // Выход из цикла
					}
				} else {

					shift_freq -= (int)(slope_real_2/0.02*1000);
				}
			}
	//        printf("shift_freq = %d ", shift_freq);
		}

		if (count_break<15){
	//        printf("shift_freq = %d ", shift_freq);
			// Вычислене фазы в серединах символов
			int len_count_phase = 0;
			for (int i = 0; i < times_change_len - 1; i++) {
				int index = times_change_phase[i] + times_change_phase_diff[i];
				count_phase[len_count_phase++] = phase_base[index];
			}

			// Вычисление разностей count_phase_diff (аналог np.diff)
			for (int i = 0; i < len_count_phase - 1; i++) {
				count_phase_diff[i] = count_phase[i + 1] - count_phase[i];
			}
			// Замена отрицательных значений в count_phase_diff на положительные
			for (int i = 0; i < len_count_phase - 1; i++) {
				if (count_phase_diff[i] < 0) {
					count_phase_diff[i] *= -1;  // Умножение на -1
				}
			}

			if (len_count_phase>200)
				len_count_phase=200;

			sort_array(count_phase_diff, len_count_phase-1);


			Point points[201];
			for (int i = 0; i < len_count_phase-1; i++) {
				points[i].value = count_phase_diff[i];
				points[i].label = 0; // Все точки изначально не классифицированы
			}
			double eps = 3;
			int min_samples = 5;
			dbscan(points, len_count_phase-1, eps, min_samples);
	//        for (int i = 0; i < 100; i++) {
	//            printf("Point %d: Cluster %d, Value=%f \n", i, points[i].label,points[i].value);
	//        }

			int max_label = find_max_label(points, len_count_phase-2);

	//        printf("Maximum label: %d\n", max_label);
			if (max_label<=5 && max_label>2){
				xil_printf("QPSK MODULATIONS\n");
			}
			else if (max_label<=2 ){
				xil_printf("BPSK MODULATIONS\n");
			}
			else if (max_label<=9 && max_label>5){

				xil_printf("8PSK MODULATIONS\n");
			}
			else xil_printf("unknown MODULATIONS\n");
		}
	//		xil_printf("ENDDDD MODULATIONS!");
	}
	else  xil_printf("NOT PSK MODULATIONS!\n");

}

