#ifndef USR_STATS_WRAPPER_H
#define USR_STATS_WRAPPER_H
// function prototypes/external variables
#ifdef __cplusplus
extern "C" {
#endif
       void open_stats_(void);
       void close_stats_(void);
       void compute_stats_(void);
       void write_restart_stats_(void);
       void read_restart_stats_(void);
#ifdef __cplusplus
}
#endif

// function wrappers

inline void mpt_open_stats(){open_stats_();}
inline void mpt_close_stats(){close_stats_();}
inline void mpt_compute_stats(){compute_stats_();}
inline void mpt_write_restart_stats(){write_restart_stats_();}
inline void mpt_read_restart_stats(){read_restart_stats_();}

#endif //USR_STATS_WRAPPER_H
