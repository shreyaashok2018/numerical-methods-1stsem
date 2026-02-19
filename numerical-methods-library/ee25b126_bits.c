/**
 * @file ee25b126_bits.c
 * @brief Implementation of Bit Operations and Hamming Distance Search
 * @author EE25B126
 * @date November 2025
 */

#include "ee25b126_ee1103.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h> 

/* ============================================================================
   9. HAMMING DISTANCE & BIT OPERATIONS
   ============================================================================ */

int hamming_distance(const unsigned char *bits1, const unsigned char *bits2, int num_bits) {
    if (!bits1 || !bits2 || num_bits <= 0) return 0;

    int dist = 0;
    int full_bytes = num_bits / 8;
    int rem_bits = num_bits % 8;

    for (int i = 0; i < full_bytes; i++) {
        unsigned char xor_byte = bits1[i] ^ bits2[i];
        dist += __builtin_popcount(xor_byte);
    }

    if (rem_bits > 0) {
        // MSB mask:
        unsigned char mask = ((1u << rem_bits) - 1) << (8 - rem_bits);
        // e.g., rem_bits=3: ((1<<3)-1) << (8-3) = (7) << 5 = 11100000

        unsigned char xor_rem = (bits1[full_bytes] ^ bits2[full_bytes]) & mask;
        dist += __builtin_popcount(xor_rem);
    }

    return dist;
}

int get_bit(const unsigned char *arr, int bit_index) {
    if (!arr || bit_index < 0) return -1;
    int byte_idx = bit_index / 8;
    int bit_in_byte = bit_index % 8;
    
    // Check if bit is set (MSB-to-LSB, e.g., 10000000 is bit 0)
    if (arr[byte_idx] & (1 << (7 - bit_in_byte))) {
        return 1;
    } else {
        return 0;
    }
}

void set_bit(unsigned char *arr, int bit_index, int value) {
    if (!arr || bit_index < 0) return;
    int byte_idx = bit_index / 8;
    int bit_in_byte = bit_index % 8;
    
    if (value == 1) {
        // Set bit to 1
        arr[byte_idx] |= (1 << (7 - bit_in_byte));
    } else {
        // Set bit to 0
        arr[byte_idx] &= ~(1 << (7 - bit_in_byte));
    }
}

void flip_bit(unsigned char *arr, int bit_index) {
    if (!arr || bit_index < 0) return;
    int byte_idx = bit_index / 8;
    int bit_in_byte = bit_index % 8;
    
    // XOR with the bit mask to flip it
    arr[byte_idx] ^= (1 << (7 - bit_in_byte));
}

void generate_random_bits(unsigned char *arr, int num_bits) {
    if (!arr || num_bits <= 0) return;
    int num_bytes = (num_bits + 7) / 8;
    memset(arr, 0, num_bytes);
    
    for (int i = 0; i < num_bits; i++) {
        if (rand() % 2 == 1) {
            set_bit(arr, i, 1);
        }
    }
}

void flip_bits_with_probability(unsigned char *arr, int num_bits, double probability) {
    if (!arr || num_bits <= 0) return;
    
    for (int i = 0; i < num_bits; i++) {
        double r = (double)rand() / RAND_MAX;
        if (r < probability) {
            flip_bit(arr, i);
        }
    }
}

void string_to_bits(const char *str, unsigned char *arr, int length) {
    if (!str || !arr || length <= 0) return;
    
    int num_bytes = (length + 7) / 8;
    memset(arr, 0, num_bytes);
    
    for (int i = 0; i < length; i++) {
        if (str[i] == '1') {
            set_bit(arr, i, 1);
        }
        // No 'else' needed, array is already zeroed
    }
}

void write_bits_to_file(const char *filename, const unsigned char *arr, int num_bits) {
    if (!filename || !arr || num_bits <= 0) return;
    
    FILE *fp = fopen(filename, "wb");
    if (!fp) return;
    
    // Write all bytes, including the last partial one 
    int num_bytes_to_write = (num_bits + 7) / 8;
    
    fwrite(arr, 1, num_bytes_to_write, fp);
    fclose(fp);
}

void read_bits_from_file(const char *filename, unsigned char *arr, int num_bits) {
    if (!filename || !arr || num_bits <= 0) return;
    
    FILE *fp = fopen(filename, "rb");
    if (!fp) return;
    
    // Read all bytes, including the last partial one 
    int num_bytes_to_read = (num_bits + 7) / 8;
    
    // Zero the buffer first
    memset(arr, 0, num_bytes_to_read);
    
    fread(arr, 1, num_bytes_to_read, fp);
    fclose(fp);
    
    // Clean trailing bits in the last byte if they matter
    int rem_bits = num_bits % 8;
    if (rem_bits > 0) {
        unsigned char mask = ((1u << rem_bits) - 1) << (8 - rem_bits);
        arr[num_bytes_to_read - 1] &= mask;
    }
}

int find_codeword(const unsigned char *bitstream, int stream_length,
                  const unsigned char *codeword, int codeword_length,
                  HammingResult *results, int max_results) {
    
    if (!bitstream || !codeword || !results || stream_length < codeword_length || max_results <= 0) {
        return 0;
    }

    int search_space = stream_length - codeword_length + 1;
    int code_bytes = (codeword_length + 7) / 8;

    // Initialize results
    for (int i = 0; i < max_results; i++) {
        results[i].location = -1;
        results[i].distance = INT_MAX;
    }

    unsigned char *window = malloc(code_bytes);
    if (!window) return 0;

    for (int pos = 0; pos < search_space; pos++) {
        // Extract window
        memset(window, 0, code_bytes);
        for (int i = 0; i < codeword_length; i++) {
            int bit = get_bit(bitstream, pos + i);
            if (bit >= 0) {
                set_bit(window, i, bit);
            }
        }

        int dist = hamming_distance(codeword, window, codeword_length);

        // Insert into top-K list (this is a simple, unoptimized insertion sort)
        if (dist < results[max_results - 1].distance) {
            results[max_results - 1].location = pos;
            results[max_results - 1].distance = dist;

            // Bubble the new result up to its correct spot
            for (int j = max_results - 1; j > 0; j--) {
                if (results[j].distance < results[j - 1].distance) {
                    HammingResult temp = results[j];
                    results[j] = results[j - 1];
                    results[j - 1] = temp;
                } else {
                    break; // In correct place
                }
            }
        }
    }

    free(window);
    
    // Count how many valid results we found
    int count = 0;
    for(int i=0; i<max_results; i++) {
        if (results[i].distance != INT_MAX) {
            count++;
        }
    }
    return count;
}

