cimport cython

from cpython.buffer cimport \
    PyBUF_SIMPLE, PyBUF_WRITABLE, \
    PyObject_CheckBuffer, PyObject_GetBuffer, PyBuffer_Release

import string
import itertools

complement = string.maketrans('ATCGN', 'TAGCN')

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

cdef inline int dont_match(char a, char b):
    if a == b: return 0
    return 1

cdef inline int dont_match_RC(char a, char b):
    if a == 'A' and b == 'T': return 0
    elif a == 'C' and b == 'G': return 0
    elif a == 'G' and b == 'C': return 0
    elif a == 'T' and b == 'A': return 0
    elif a == 'N' and b == 'N': return 0
    else: return 1
    

# Align with mismatch, find first and move on, assumes only one
#@cython.boundscheck(False)
def fuzz_align(str left_seq not None, 
               str right_seq not None, 
               str adapter not None, 
               int mismatch, 
               int MIN_MATCH_LENGTH):
    """Find the optimal alignment of left_seq and right seq, 
    where the right portion of left seq is allowed to overlap
    the left portion of right seq. 

    """
    cdef int curr_offset = -1
    cdef float curr_mismatch = mismatch

    cdef int left_seq_len = len(left_seq)    
    cdef Py_buffer left_seq_view
    PyObject_GetBuffer(left_seq.upper(), &left_seq_view, PyBUF_SIMPLE)
    cdef  char* c_left_seq = < char *>left_seq_view.buf
    
    cdef int right_seq_len = len(right_seq)    
    cdef Py_buffer right_seq_view
    PyObject_GetBuffer(right_seq.upper(), &right_seq_view, PyBUF_SIMPLE)
    cdef  char* c_right_seq = < char *>right_seq_view.buf

    cdef int adapter_len = len(adapter)    
    cdef Py_buffer adapter_view
    PyObject_GetBuffer(adapter.upper(), &adapter_view, PyBUF_SIMPLE)
    cdef  char* c_adapter = < char *>adapter_view.buf
    
    #left_seq = left_seq.upper()
    #right_seq = right_seq.upper()
    #adapter = adapter.upper()

    cdef float dist
    cdef int offset
    cdef int i
    # loop through all allowable offsets
    #for offset in xrange(0, int_max(1, left_seq_len - MIN_MATCH_LENGTH)): 
    for offset in xrange(left_seq_len - MIN_MATCH_LENGTH): #xrange(adapter_len+1):
        dist = 0
        # if the reads overlap, then the reverse complemenet of the left 
        # side of the left sequence should match the left side of the right 
        # sequence. Calculate the edit distance
        for i in xrange(left_seq_len-offset):
            if dont_match_RC(c_left_seq[left_seq_len-offset-1-i], 
                             c_right_seq[i]):
                dist += 1
        # In addition, if the reads overlap, then the adapters sequence should
        # overlap both reads. 
        if offset <= adapter_len:
            for i in xrange(offset):
                #if left_seq[offset-i-1].translate(complement) != adapter[i]:
                if dont_match(c_left_seq[left_seq_len-offset+i], c_adapter[i]):
                    dist += 0.5
                if dont_match(c_right_seq[right_seq_len-offset+i], c_adapter[i]):
                    dist += 0.5
        
        if dist <= curr_mismatch:
            curr_offset = offset
            curr_mismatch = dist

    """
    print left_seq
    print right_seq
    
    print left_seq[-curr_offset:]
    print right_seq[-curr_offset:]

    print left_seq[:-curr_offset]
    print right_seq[:-curr_offset]
    """

    PyBuffer_Release(&left_seq_view)
    PyBuffer_Release(&right_seq_view)
    PyBuffer_Release(&adapter_view)
    if curr_offset == -1: return None
    return curr_offset, curr_mismatch
