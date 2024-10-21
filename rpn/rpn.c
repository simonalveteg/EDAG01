#include <stdio.h>
#include <ctype.h>
#define N		     (10)
#define ADD      '+'
#define SUB      '-'
#define MUL      '*'
#define DIV      '/'

int is_binary_operator(int c) {
  return c == ADD || c == SUB || c == MUL || c == DIV;
}

int stack_size(int top_of_stack) {
  return top_of_stack + 1;
}

void print_error(int line, char error) {
  // Make sure \n gets printed correctly by converting char to string.
  char *error_string = (error == '\n') ? "\\n" : (char[2]){error, '\0'};
  printf("line %d: error at %s\n", line, error_string);
}

int main(void)
{
  int stack[N];
  int c;
  int line = 0;
  int top_of_stack = -1; // keep track of top of stack.
  int previous_was_digit = 0;
  int error_found = '\0';

  while ((c = getchar()) != EOF) {
    if (error_found != '\0') {
      if (c == '\n') {
        line++;
        print_error(line, (char)error_found);
        error_found = '\0';
        top_of_stack = -1; // reset stack.
      }
      
      continue;
    }
    if (isdigit(c)) {
      int d = c - '0'; // Convert char to int.
      // If previous char was digit, this is still the same number.
      if (previous_was_digit) {
        stack[top_of_stack] = 10*stack[top_of_stack] + d;
      } else {
        // Previous wasn't digit, this is a new number.
        // Error if stack is full.
        if (stack_size(top_of_stack) >= N) {
          error_found = c; 
          continue;
        }
        top_of_stack++;
        stack[top_of_stack] = d;
      }
      previous_was_digit = 1;
    } else if (is_binary_operator(c)) {
      previous_was_digit = 0;
      if (stack_size(top_of_stack) < 2) {
        error_found = c;
        continue;
      }
      int a = stack[top_of_stack];
      top_of_stack--;
      int b = stack[top_of_stack];
      switch (c) {
        case ADD:
          stack[top_of_stack] = b + a;
          break;
        case SUB:
          stack[top_of_stack] = b - a;
          break;
        case MUL:
          stack[top_of_stack] = b * a;
          break;
        case DIV:
          if (a == 0) {
            error_found = DIV;
            continue;
          }
          stack[top_of_stack] = b / a;
          break;
      }
    } else {
      previous_was_digit = 0;
      switch (c) {
        case ' ': break; // do nothing.
        case '\n': 
          // End of line reached, print result.
          line++;
          if (stack_size(top_of_stack) != 1) {
            print_error(line, (char)c);
            top_of_stack = -1; // reset stack.
            break;
          }
          printf("line %d: %d\n", line, stack[top_of_stack]);
          top_of_stack--;
          break;
        default:
          // Invalid input, throw error.
          error_found = c;
      }  
    }
  }

	return 0;
}
