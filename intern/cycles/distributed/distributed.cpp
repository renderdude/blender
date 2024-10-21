#include "distributed.h"
#include <sstream>

#include "mpl/comm_group.hpp"

// And ZeroMQ
#include "zmq.hpp"

CCL_NAMESPACE_BEGIN

Distributed::Distributed(std::string bind_to, std::string connect_to)
    : comm_world{mpl::environment::comm_world()}, is_distributed(false)
{
  char hostname[256];
  if (gethostname(hostname, 256)) {
    printf("Unable to get hostname\n");
    memcpy(&hostname, "UNKNOWN", 7);
  }

  if (!connect_to.empty() || !bind_to.empty()) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib;

    std::stringstream host_string;
    // Set the hostname to the host plus a random number to distinguish
    // amongst processes on the same host
    host_string << hostname << ":" << distrib(gen);

    zmq::context_t context(1);

    // get a reference to communicator "world"

    int size = comm_world.size();
    int rank = comm_world.rank();

    mpl::inter_communicator *inter_comm;

    if (rank == 0) {
      if (!connect_to.empty()) {
        zmq::socket_t socket(context, zmq::socket_type::req);
        socket.connect(connect_to);
        // Let the server know that we're here
        zmq::message_t msg(host_string.str());
        socket.send(msg, zmq::send_flags::none);
        // Now, get the MPI port information for connecting
        zmq::message_t recv_msg;
        zmq::recv_result_t result = socket.recv(recv_msg);
        std::string port_info(static_cast<char *>(recv_msg.data()), recv_msg.size());
        std::cout << host_string.str() << ": connecting" << std::endl;
        inter_comm = new mpl::inter_communicator(comm_world.connect(port_info, mpl::info(), 0));
        socket.close();
        // Setup the new intercomm with the 'client' having the high-process
        // numbers
        inter_comm->barrier();
        inter_comm_world = mpl::communicator(*inter_comm, mpl::communicator::order_high);
      }
      else {
        zmq::socket_t socket(context, zmq::socket_type::rep);
        socket.bind(bind_to);
        mpl::port port;
        if (port.open() != MPI_SUCCESS) {
          std::cerr << "open failed\n";
          exit(1);
        }
        // Wait for the client to connect
        zmq::message_t recv_msg;
        zmq::recv_result_t result = socket.recv(recv_msg);
        // Now send the MPI port info
        zmq::message_t msg(port.name());
        socket.send(msg, zmq::send_flags::none);
        std::cerr << host_string.str() << ": Waiting for connection." << std::endl;
        inter_comm = new mpl::inter_communicator(comm_world.accept(port.name(), mpl::info(), 0));
        std::cerr << "Client connected." << std::endl;
        socket.close();
        // The server has rank 0 in the inter-comm
        inter_comm->barrier();
        inter_comm_world = mpl::communicator(*inter_comm, mpl::communicator::order_low);
      }

      std::cerr << "intra communicator has size " << inter_comm_world.size()
                << " with rank = " << inter_comm_world.rank() << std::endl;
    }

    // MPI takes over from here on out
    inter_comm_world.barrier();
  }
}

CCL_NAMESPACE_END
