#include "distributed.h"
#include <iostream>
#include <random>
#include <sstream>

#include "mpl/comm_group.hpp"

CCL_NAMESPACE_BEGIN

Distributed::Distributed(bool reverse_connect)
    : comm_world{mpl::environment::comm_world()}, is_render_server(reverse_connect)

{
  char hostname[256];
  if (gethostname(hostname, 256)) {
    printf("Unable to get hostname\n");
    strcpy(hostname, "UNKNOWN");
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib;

  std::stringstream host_string;
  // Set the hostname to the host plus a random number to distinguish
  // amongst processes on the same host
  host_string << hostname << ":" << distrib(gen);

  int size = comm_world.size();
  int rank = comm_world.rank();

  if (rank == 0) {
    if (is_render_server) {
      mpl::port port;
      port.lookup_name("papillon");
      std::cout << host_string.str() << ": connecting" << std::endl;
      _inter_comm = new mpl::inter_communicator(comm_world.connect(port.name(), mpl::info(), 0));
      // Setup the new intercomm with the 'client' having the high-process numbers
      _inter_comm->barrier();
      inter_comm_world = mpl::communicator(*_inter_comm, mpl::communicator::order_high);
    }
    else {  // Server side
      mpl::info info;
      info.set("ompi_global_scope", "true");

      mpl::port port(info);
      if (port.open() != MPI_SUCCESS) {
        std::cerr << "open failed\n";
        exit(1);
      }
      port.publish_name("papillon");
      // Wait for the client to connect
      std::cerr << host_string.str() << ": Waiting for connection." << std::endl;
      _inter_comm = new mpl::inter_communicator(comm_world.accept(port.name(), mpl::info(), 0));
      std::cerr << "Client connected." << std::endl;
      // The server has rank 0 in the inter-comm
      _inter_comm->barrier();
      inter_comm_world = mpl::communicator(*_inter_comm, mpl::communicator::order_low);
      port.unpublish_name("papillon");

    }

    // MPI takes over from here on out
    inter_comm_world.barrier();
  }
}

CCL_NAMESPACE_END
