cdef extern from "infernal.h" nogil:

    void  init_ilogsum()
    int   ILogsum(int s1, int s2)
    int   ILogsumNI(int s1, int s2)
    int   ILogsumNI_diff(int s1a, int s1b, int s2a, int s2b, int db)
    void  FLogsumInit()
    float LogSum2(float p1, float p2)
    float FLogsum(float p1, float p2)