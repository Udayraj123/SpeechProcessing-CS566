// MUST CONFIGURE THESE BEFORE RUNNING ON NEW DEVICE
#define DC_SHIFT 		 0
#define CURR_DIR "C:/Users/Udayraj021/Downloads/coding/SpeechProcessing-CS566/Digits/"
#define INPUT_TRIMMED_ALREADY 	 1

// User Config
#define SHOW_RIs 	 0
#define SHOW_AIs 	 0
#define SHOW_CIs 	 1
#define SAVE_TRIMMED_CLIPS 	 0

// Model Parameters 
#define P_ORDER 		 12
#define USE_HAMMING 	 1
#define USE_TOKHURA	 	 1
#define WINDOW_SIZE 	 320
#define WINDOW_STRIDE	 80
// #define FIRST_N_FRAMES 	 5  <- same as OBSERVATIONS_LIM

// Default Config
#define DEFAULT_DURN	"5"
#define N_AMP 			5000 
#define N_CLIPS_LIMIT 	15 // avoid too many file writes
#define SAMPLERATE 		16000.0

#define M_SILENCE_ENERGY 	22500
#define SILENCE_WINDOW 	 3290

// LBG
// Turn this on if C0 is present
#define SKIP_FIRST_CI 1
#define CODEBOOK_SIZE 	 16
#define KMEANS_ITERATIONS 	 12
// EM Algorithm
#define OBSERVATIONS_LIM 		 100
#define NUM_ITERATIONS 		 	 20
// No of epochs to average A,B matrices and rerun EM algo
#define NUM_EPOCHS 3