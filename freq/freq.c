#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#define WORD_SIZE 5

typedef struct word_t {
    int count;
    char word[WORD_SIZE];
    struct word_t* next;
} word_t;

typedef struct list_t {
    word_t* first;
    word_t* last;
} list_t;

bool is_prime(int n) {
  if (n <= 1) return false;         // 0 and 1 are not prime numbers
  if (n <= 3) return true;          // 2 and 3 are prime numbers
  if (n % 2 == 0 || n % 3 == 0) return false;  // Exclude multiples of 2 and 3

  for (int i = 5; i * i <= n; i += 6) {
      if (n % i == 0 || n % (i + 2) == 0) return false;
  }

  return true;
}

void free_list(list_t* list) {
  if (list == NULL) return;
  word_t *word = list->first;
  word_t *next;

  while (word != NULL) {
    next = word->next;
    word->next = NULL;
    free(word);
    word = next;
  }
  list->first = list->last = NULL;
  free(list);
}

void add_word(list_t* list, char* word) {
  // Create a new word and add it to the list.
  word_t* new_word = malloc(sizeof(word_t));
  new_word->count = 1;
  strcpy(new_word->word, word);
  new_word->next = NULL;

  // Add first if no words are present.
  if (list->first == NULL && list->last == NULL) {
    list->first = new_word;
    list->last = new_word;
    printf("added %s\n", word);
    return;
  } 

  // update word with same word if possible.
  word_t *current = list->first;
  while (current != NULL) {
    if (strcmp(current->word, word) == 0) {
      current->count++;
      free(new_word);
      printf("counted %s\n", word);
      return;
    }
    current = current->next;
  }
  
  // Add word last.
  list->last->next = new_word;
  list->last = new_word;

  printf("added %s\n", word);
}

void delete_word(list_t* list, char* word) {
  word_t *current = list->first;
  word_t *previous = NULL;
  printf("trying to delete %s: ", word);
  while (current != NULL) {
    if (strcmp(current->word, word) == 0) {
      // Word match found.
      if (list->first == list->last) {
        // Word was only item.
        list->last = list->first = NULL;
      } else if (current == list->first) {
        // Word found in first item.
        list->first = current->next;
      } else if (current == list->last) {
        // Word found in last item.
        previous->next = NULL;
        list->last = previous;
      } else {
        previous->next = current->next;
      }
      printf("deleted\n");
      free(current);
      return;
    }
    previous = current;
    current = current->next;
  }
  printf("not found\n");
}

int main() 
{
  char str[WORD_SIZE];
  int line = 0;

  list_t *list = malloc(sizeof(list_t));
  list->first = list->last = NULL;

  while(scanf("%s", str) != EOF) {
    line++;
    //printf("Word read: %s, line: %d, prime: %d\n", str, line, is_prime(line));
    // Remove if line is prime, count otherwise.
    if (is_prime(line)) {
      delete_word(list, str);
    } else {
      add_word(list, str);
    }
  }
  
  word_t *current = list->first;
  word_t *max = NULL;
  while (current != NULL) {
    if (max == NULL || current->count > max->count) {
      max = current;
    }
    current = current->next;
  }
  printf("result: %s %d\n", max->word, max->count);
  free_list(list);

  return 0;
}
