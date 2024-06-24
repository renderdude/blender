#ifndef DISTRIBUTED_H
#define DISTRIBUTED_H

#include <random>
#include <vector>

#include "mpl/mpl.hpp"
#include "mpl/comm_group.hpp"

// And ZeroMQ
#include "zmq.hpp"

#include "app/cycles_standalone.h"

CCL_NAMESPACE_BEGIN

/** @brief distributed.
 * @details
 */

class Distributed {
public:
  /// @name Initialization
  ///@{
  Distributed(Options &options)
      : comm_world{mpl::environment::comm_world()}, is_distributed(false) {
    char hostname[256];
    if (gethostname(hostname, 256)) {
      printf("Unable to get hostname\n");
      memcpy(&hostname, "UNKNOWN", 7);
    }

    if (!options.connect_to.empty() || !options.bind_to.empty()) {
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
        if (!options.connect_to.empty()) {
          zmq::socket_t socket(context, zmq::socket_type::req);
          socket.connect(options.connect_to);
          // Let the server know that we're here
          zmq::message_t msg(host_string.str());
          socket.send(msg, zmq::send_flags::none);
          // Now, get the MPI port information for connecting
          zmq::message_t recv_msg;
          zmq::recv_result_t result = socket.recv(recv_msg);
          std::string port_info(static_cast<char *>(recv_msg.data()),
                                recv_msg.size());
          std::cout << host_string.str() << ": connecting" << std::endl;
          inter_comm = new mpl::inter_communicator(
              comm_world.connect(port_info, mpl::info(), 0));
          socket.close();
          // Setup the new intercomm with the 'client' having the high-process
          // numbers
          inter_comm->barrier();
          inter_comm_world =
              mpl::communicator(*inter_comm, mpl::communicator::order_high);
        } else {
          zmq::socket_t socket(context, zmq::socket_type::rep);
          socket.bind(options.bind_to);
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
          std::cerr << host_string.str() << ": Waiting for connection."
                    << std::endl;
          inter_comm = new mpl::inter_communicator(
              comm_world.accept(port.name(), mpl::info(), 0));
          std::cerr << "Client connected." << std::endl;
          socket.close();
          // The server has rank 0 in the inter-comm
          inter_comm->barrier();
          inter_comm_world =
              mpl::communicator(*inter_comm, mpl::communicator::order_low);
        }

        std::cerr << "intra communicator has size " << inter_comm_world.size()
                  << " with rank = " << inter_comm_world.rank() << std::endl;
      }

      // MPI takes over from here on out
      inter_comm_world.barrier();
    }
  }
  ///@}

  /// @name Duplication
  ///@{
  Distributed(Distributed const &) = default;
  Distributed &operator=(Distributed const &) = default;
  ///@}

  /// @name Move
  ///@{
  Distributed(Distributed &&) = default;
  Distributed &operator=(Distributed &&) = default;
  ///@}

  /// @name Destruction
  ///@{
  ~Distributed() = default;
  ///@}

  /// @name Access
  ///@{
  mpl::communicator comm_world, inter_comm_world;
  ///@}
  /// @name Measurement
  ///@{
  ///@}
  /// @name Comparison
  ///@{
  ///@}
  /// @name Status report
  ///@{
  bool is_distributed;
  ///@}
  /// @name Status setting
  ///@{
  ///@}
  /// @name Cursor movement
  ///@{
  ///@}
  /// @name Element change
  ///@{
  ///@}
  /// @name Removal
  ///@{
  ///@}
  /// @name Resizing
  ///@{
  ///@}
  /// @name Transformation
  ///@{
  ///@}
  /// @name Conversion
  ///@{
  std::string to_string() const;
  ///@}
  /// @name Basic operations
  ///@{
  ///@}
  /// @name Miscellaneous
  ///@{
  ///@}
  /// @name Obsolete
  ///@{
  ///@}
  /// @name Inapplicable
  ///@{
  ///@}

protected:
  ///@{
  ///@}

private:
  /* data */

}; // end of class Distributed
CCL_NAMESPACE_END

#endif // DISTRIBUTED_H
