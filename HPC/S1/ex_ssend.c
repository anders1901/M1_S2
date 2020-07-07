#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>


#define SIZE_H_N 50

int main(int argc, char *argv[]) {
  int my_rank, p, source, dest;
  int tag = 0;
  char message[100];
  MPI_Status status;
  char hostname[SIZE_H_N];

  gethostname(hostname, SIZE_H_N);

  /* Initialisation */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,   &p);
  printf("\nAVANT BOUCLE IF \n - my_rank = %d\n - p = %d\n", my_rank, p);
  if (my_rank != 0) {
    /* Création du message */
    printf("BOUCLE IF \n - my_rank = %d\n", my_rank);
    sprintf(message, "Coucou du processus #%d depuis %s !", my_rank, hostname);
    dest = my_rank + 1;

    MPI_Send(message, strlen(message) + 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
  }else{
    printf("BOUCLE ELSE \n - my_rank = %d\n", my_rank);
    for (source = 1; source < p; source++) {
      printf("BOUCLE FOR \n - my_rank = %d\n - source = %d\n", my_rank, source);
      MPI_Recv(message, 100, MPI_CHAR, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
      printf("Sur %s, le processus #%d a recu le message : %s\n", hostname, my_rank, message);
    }
  }

  /* Desactivation */
  MPI_Finalize();
  return 0;
}
