// MUST CONFIGURE THESE BEFORE RUNNING ON NEW DEVICE
#define DC_SHIFT 		 0
#define CURR_DIR "C:/Users/Udayraj021/Downloads/coding/SpeechProcessing-CS566/Vowels/"
#define INPUT_FOLDER "input_data/"
#define INPUT_TRIMMED_ALREADY 	 0
#define TEST_TRIMMED_ALREADY 	 1
#define DEFAULT_RECORD_FILE "record.txt"

// User Config
#define SHOW_RIs 	 0
#define SHOW_AIs 	 0
#define SHOW_CIs 	 1
#define SAVE_TRIMMED_CLIPS 	 1

// Model Parameters 
#define P_ORDER 		 12
#define USE_HAMMING 	 1
#define USE_TOKHURA	 	 1
#define WINDOW_SIZE 	 320
#define WINDOW_STRIDE	 320
#define FIRST_N_FRAMES 	 5 

// Default Config
#define DEFAULT_DURN	"5"
#define N_AMP 			5000 
#define N_CLIPS_LIMIT 	15 // avoid too many file writes
#define SAMPLERATE 		16000.0

#define M_SILENCE_ENERGY 	22500
#define SILENCE_WINDOW 	 3790