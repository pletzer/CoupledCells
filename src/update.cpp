#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "macros.h"
#include "computelib.h"

using namespace std;
extern time_stamps t_stamp;

void determine_source_destination(grid_parms grid, int source[], int dest[])
{
	if (grid.nbrs[local][UP] >= 0)
	{
		dest[UP] = grid.nbrs[local][UP];
		source[UP] = grid.nbrs[local][UP];
	}
	else if (grid.nbrs[local][UP] < 0)
	{
		dest[UP] = grid.rank_branch; //MPI_PROC_NULL;
		source[UP] = grid.rank_branch; //MPI_PROC_NULL;
	}

	if (grid.nbrs[local][DOWN] >= 0)
	{
		dest[DOWN] = grid.nbrs[local][DOWN];
		source[DOWN] = grid.nbrs[local][DOWN];
	}
	else if (grid.nbrs[local][DOWN] < 0)
	{
		dest[DOWN] = grid.rank_branch; //MPI_PROC_NULL;
		source[DOWN] = grid.rank_branch; //MPI_PROC_NULL;
	}

	if (grid.nbrs[local][LEFT] >= 0)
	{
		dest[LEFT] = grid.nbrs[local][LEFT];
		source[LEFT] = grid.nbrs[local][LEFT];
	}
	else if (grid.nbrs[local][LEFT] < 0)
	{
		dest[LEFT] = grid.rank_branch; //MPI_PROC_NULL;
		source[LEFT] = grid.rank_branch; //MPI_PROC_NULL;
	}

	if (grid.nbrs[local][RIGHT] >= 0)
	{
		dest[RIGHT] = grid.nbrs[local][RIGHT];
		source[RIGHT] = grid.nbrs[local][RIGHT];
	}
	else if (grid.nbrs[local][RIGHT] < 0)
	{
		dest[RIGHT] = grid.rank_branch; //MPI_PROC_NULL;
		source[RIGHT] = grid.rank_branch; // MPI_PROC_NULL;
	}
}

grid_parms communicate_num_recv_elements_to_nbrs(grid_parms grid)
{
	int err;
	MPI_Request reqs[8];
	MPI_Status stats[8];

	int source[4], dest[4];
	int tag;
	int num[4];
	determine_source_destination(grid, source, dest);

	tag = 0;
	int tmp[4];
	CHECK_MPI_ERROR(MPI_Irecv(&grid.num_elements_recv_up, 1, MPI_INT, source[UP], tag, grid.cart_comm, &reqs[4 + UP]));
	CHECK_MPI_ERROR(MPI_Isend(&grid.num_elements_send_up, 1, MPI_INT, dest[UP], tag, grid.cart_comm, &reqs[UP]));
	CHECK_MPI_ERROR(MPI_Irecv(&grid.num_elements_recv_down, 1, MPI_INT, source[DOWN], tag, grid.cart_comm, &reqs[4 + DOWN]));
	CHECK_MPI_ERROR(MPI_Isend(&grid.num_elements_send_down, 1, MPI_INT, dest[DOWN], tag, grid.cart_comm, &reqs[DOWN]));
	CHECK_MPI_ERROR(MPI_Irecv(&grid.num_elements_recv_left, 1, MPI_INT, source[LEFT], tag, grid.cart_comm, &reqs[4 + LEFT]));
	CHECK_MPI_ERROR(MPI_Isend(&grid.num_elements_send_left, 1, MPI_INT, dest[LEFT], tag, grid.cart_comm, &reqs[LEFT]));
	CHECK_MPI_ERROR(MPI_Irecv(&grid.num_elements_recv_right, 1, MPI_INT, source[RIGHT], tag, grid.cart_comm, &reqs[4 + RIGHT]));
	CHECK_MPI_ERROR(MPI_Isend(&grid.num_elements_send_right, 1, MPI_INT, dest[RIGHT], tag, grid.cart_comm, &reqs[RIGHT]));

	MPI_Waitall(8, reqs, stats);

	return (grid);
}

// Carry out asynchronous communication send and receive of edge cell data.
void communication_async_send_recv(grid_parms grid, double** sendbuf, double** recvbuf, SMC_cell** smc, EC_cell** ec)
{
	/// 8 MPI variables of types Request and Status are declared.
	/// 4 of each are for sending information to 4 neighbours and the other 4 are to retrive the receive operation status.
	MPI_Request reqs[8];
	MPI_Status stats[8];

	/// Two arrays, source and dest, hold the ranks for the communicating tasks.
	/// dest has ranks to which my task will send a message to.
	/// source contains tasks from which I expect a message.
	int source[4], dest[4];

	// Get nearest neighbours indices.
	determine_source_destination(grid, source, dest);

	t_stamp.update_sendbuf_t1 = MPI_Wtime();

	// Prepare the buffer for exchanging edge cell data with ghost cells.
	communication_update_sendbuf(grid, sendbuf, smc, ec);

	t_stamp.update_sendbuf_t2 = MPI_Wtime();
	t_stamp.diff_update_sendbuf = t_stamp.update_sendbuf_t2 - t_stamp.update_sendbuf_t1;

	t_stamp.async_comm_calls_t1 = MPI_Wtime();

	/// Communication block
	CHECK_MPI_ERROR(MPI_Irecv(&recvbuf[UP][0], grid.num_elements_recv_up, MPI_DOUBLE, source[UP], MPI_ANY_TAG, grid.cart_comm, &reqs[4 + UP]));

	CHECK_MPI_ERROR(MPI_Irecv(&recvbuf[DOWN][0], grid.num_elements_recv_down, MPI_DOUBLE, source[DOWN], MPI_ANY_TAG, grid.cart_comm, &reqs[4 + DOWN]));

	CHECK_MPI_ERROR(MPI_Irecv(&recvbuf[LEFT][0], grid.num_elements_recv_left, MPI_DOUBLE, source[LEFT], MPI_ANY_TAG, grid.cart_comm, &reqs[4 + LEFT]));

	CHECK_MPI_ERROR(MPI_Irecv(&recvbuf[RIGHT][0], grid.num_elements_recv_right, MPI_DOUBLE, source[RIGHT], MPI_ANY_TAG, grid.cart_comm, &reqs[4 + RIGHT]));

	CHECK_MPI_ERROR(MPI_Isend(sendbuf[UP], grid.num_elements_send_up, MPI_DOUBLE, dest[UP], 0, grid.cart_comm, &reqs[UP]));

	CHECK_MPI_ERROR(MPI_Isend(sendbuf[DOWN], grid.num_elements_send_down, MPI_DOUBLE, dest[DOWN], 0, grid.cart_comm, &reqs[DOWN]));

	CHECK_MPI_ERROR(MPI_Isend(sendbuf[LEFT], grid.num_elements_send_left, MPI_DOUBLE, dest[LEFT], 0, grid.cart_comm, &reqs[LEFT]));

	CHECK_MPI_ERROR(MPI_Isend(sendbuf[RIGHT], grid.num_elements_send_right, MPI_DOUBLE, dest[RIGHT], 0, grid.cart_comm, &reqs[RIGHT]));

	t_stamp.async_comm_calls_t2 = MPI_Wtime();
	t_stamp.diff_async_comm_calls = t_stamp.async_comm_calls_t2 - t_stamp.async_comm_calls_t1;

	t_stamp.async_comm_calls_wait_t1 = MPI_Wtime();
	MPI_Waitall(8, reqs, stats);
	t_stamp.async_comm_calls_wait_t2 = MPI_Wtime();
	t_stamp.diff_async_comm_calls_wait = t_stamp.async_comm_calls_wait_t2 - t_stamp.async_comm_calls_wait_t1;

	int remote_source[2], remote_dest[2];
	MPI_Request remote_req[4];
	MPI_Status remote_status[4];
	int tag_remote_1 = 3, tag_remote_2 = 4;

	if (grid.my_domain.internal_info.boundary_tag == 'I')
	{
		for (int i = 0; i < 2; i++) {
			remote_source[i] = grid.nbrs[remote][i];
			remote_dest[i] = grid.nbrs[remote][i];
		}

		t_stamp.remote_async_comm_calls_t1 = MPI_Wtime();

		CHECK_MPI_ERROR(MPI_Irecv(recvbuf[UP], grid.num_elements_recv_up, MPI_DOUBLE, remote_source[UP], tag_remote_1, grid.universe,
				&remote_req[2 + UP]));

		CHECK_MPI_ERROR(MPI_Irecv(recvbuf[DOWN], grid.num_elements_recv_down, MPI_DOUBLE, remote_source[DOWN], tag_remote_1, grid.universe,
				&remote_req[2 + DOWN]));

		CHECK_MPI_ERROR(MPI_Isend(sendbuf[UP], grid.num_elements_send_up, MPI_DOUBLE, remote_dest[UP], tag_remote_1, grid.universe,
				&remote_req[UP]));

		CHECK_MPI_ERROR(MPI_Isend(sendbuf[DOWN], grid.num_elements_send_down, MPI_DOUBLE, remote_dest[DOWN], tag_remote_1, grid.universe,
				&remote_req[DOWN]));

		t_stamp.remote_async_comm_calls_t2 = MPI_Wtime();
		t_stamp.diff_remote_async_comm_calls = t_stamp.remote_async_comm_calls_t2 - t_stamp.remote_async_comm_calls_t1;

		t_stamp.remote_async_comm_calls_wait_t1 = MPI_Wtime();
		MPI_Waitall(4, remote_req, remote_status);
		t_stamp.remote_async_comm_calls_wait_t2 = MPI_Wtime();
		t_stamp.diff_remote_async_comm_calls_wait = t_stamp.remote_async_comm_calls_wait_t2 - t_stamp.remote_async_comm_calls_wait_t1;

	}
	else if ((grid.my_domain.internal_info.boundary_tag == 'T') || (grid.my_domain.internal_info.boundary_tag == 'B'))
	{
		for (int i = 0; i < 2; i++) {
			remote_source[i] = grid.nbrs[remote][i];
			remote_dest[i] = grid.nbrs[remote][i];
		}

		t_stamp.remote_async_comm_calls_t1 = MPI_Wtime();

		CHECK_MPI_ERROR(MPI_Irecv(recvbuf[UP], grid.num_elements_recv_up, MPI_DOUBLE, remote_source[UP], tag_remote_1, grid.universe,
						&remote_req[2 + UP]));

		CHECK_MPI_ERROR(MPI_Irecv(recvbuf[DOWN], grid.num_elements_recv_down, MPI_DOUBLE, remote_source[DOWN], tag_remote_1,
						grid.universe, &remote_req[2 + DOWN]));

		CHECK_MPI_ERROR(MPI_Isend(sendbuf[UP], grid.num_elements_send_up, MPI_DOUBLE, remote_dest[UP], tag_remote_1, grid.universe,
						&remote_req[UP]));

		CHECK_MPI_ERROR(MPI_Isend(sendbuf[DOWN], grid.num_elements_send_down, MPI_DOUBLE, remote_dest[DOWN], tag_remote_1,
						grid.universe, &remote_req[DOWN]));

		t_stamp.remote_async_comm_calls_t2 = MPI_Wtime();
		t_stamp.diff_remote_async_comm_calls = t_stamp.remote_async_comm_calls_t2 - t_stamp.remote_async_comm_calls_t1;

		t_stamp.remote_async_comm_calls_wait_t1 = MPI_Wtime();
		MPI_Waitall(4, remote_req, remote_status);

		t_stamp.remote_async_comm_calls_wait_t2 = MPI_Wtime();
		t_stamp.diff_remote_async_comm_calls_wait = t_stamp.remote_async_comm_calls_wait_t2 - t_stamp.remote_async_comm_calls_wait_t1;
	}

	t_stamp.update_recvbuf_t1 = MPI_Wtime();

	// Unpack received data into ghost cells.
	communication_update_recvbuf(grid, recvbuf, smc, ec);

	t_stamp.update_recvbuf_t2 = MPI_Wtime();
	t_stamp.diff_update_recvbuf = t_stamp.update_recvbuf_t2 - t_stamp.update_recvbuf_t1;
}

// Prepare the sendbuf pub putting ECs's and SMCs's edge data to be sent to the neighbours's ghost cells.
// Specific to KNBGR.
void communication_update_sendbuf(grid_parms grid, double** sendbuf, SMC_cell** smc, EC_cell** ec)
{
	int k, buf_offset;

///UP direction

	buf_offset = grid.added_info_in_send_buf;

	k = 0;

	for (int i = (int) sendbuf[UP][0]; i <= (int) sendbuf[UP][1]; i++) {
		int j = 1;
		sendbuf[UP][buf_offset + k + 0] = smc[i][j].vars[smc_Ca];
		sendbuf[UP][buf_offset + k + 1] = smc[i][j].vars[smc_Vm];
		sendbuf[UP][buf_offset + k + 2] = smc[i][j].vars[smc_IP3];
		k += grid.num_coupling_species_smc;
	}

	//Setting up the offset to transfer the EC info into SEND BUFFER
	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc * (sendbuf[UP][1] - sendbuf[UP][0] + 1);

	k = 0;
	for (int i = (int) sendbuf[UP][2]; i <= (int) sendbuf[UP][3]; i++) {
		int j = 1;
		sendbuf[UP][buf_offset + k + 0] = ec[i][j].vars[ec_Ca];
		sendbuf[UP][buf_offset + k + 1] = ec[i][j].vars[ec_Vm];
		sendbuf[UP][buf_offset + k + 2] = ec[i][j].vars[ec_IP3];
		k += grid.num_coupling_species_ec;
	}

///DOWN direction

	buf_offset = grid.added_info_in_send_buf;

	k = 0;
	for (int i = (int) sendbuf[DOWN][0]; i <= (int) sendbuf[DOWN][1]; i++) {
		int j = grid.num_smc_axially;
		sendbuf[DOWN][buf_offset + k + 0] = smc[i][j].vars[smc_Ca];
		sendbuf[DOWN][buf_offset + k + 1] = smc[i][j].vars[smc_Vm];
		sendbuf[DOWN][buf_offset + k + 2] = smc[i][j].vars[smc_IP3];
		k += grid.num_coupling_species_smc;
	}

	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc * (sendbuf[DOWN][1] - sendbuf[DOWN][0] + 1);

	k = 0;
	for (int i = (int) sendbuf[DOWN][2]; i <= (int) sendbuf[DOWN][3]; i++) {
		int j = grid.num_ec_axially;
		sendbuf[DOWN][buf_offset + k + 0] = ec[i][j].vars[ec_Ca];
		sendbuf[DOWN][buf_offset + k + 1] = ec[i][j].vars[ec_Vm];
		sendbuf[DOWN][buf_offset + k + 2] = ec[i][j].vars[ec_IP3];
		k += grid.num_coupling_species_ec;
	}

///LEFT direction

	buf_offset = grid.added_info_in_send_buf;

	k = 0;
	for (int j = (int) sendbuf[LEFT][0]; j <= (int) sendbuf[LEFT][1]; j++) {
		int i = 1;
		sendbuf[LEFT][buf_offset + k + 0] = smc[i][j].vars[smc_Ca];
		sendbuf[LEFT][buf_offset + k + 1] = smc[i][j].vars[smc_Vm];
		sendbuf[LEFT][buf_offset + k + 2] = smc[i][j].vars[smc_IP3];
		k += grid.num_coupling_species_smc;
	}

	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc * (sendbuf[LEFT][1] - sendbuf[LEFT][0] + 1);

	k = 0;
	for (int j = (int) sendbuf[LEFT][2]; j <= (int) sendbuf[LEFT][3]; j++) {
		int i = 1;
		sendbuf[LEFT][buf_offset + k + 0] = ec[i][j].vars[ec_Ca];
		sendbuf[LEFT][buf_offset + k + 1] = ec[i][j].vars[ec_Vm];
		sendbuf[LEFT][buf_offset + k + 2] = ec[i][j].vars[ec_IP3];
		k += grid.num_coupling_species_ec;
	}

///RIGHT direction

	buf_offset = grid.added_info_in_send_buf;

	k = 0;
	for (int j = (int) sendbuf[RIGHT][0]; j <= (int) sendbuf[RIGHT][1]; j++) {
		int i = grid.num_smc_circumferentially;
		sendbuf[RIGHT][buf_offset + k + 0] = smc[i][j].vars[smc_Ca];
		sendbuf[RIGHT][buf_offset + k + 1] = smc[i][j].vars[smc_Vm];
		sendbuf[RIGHT][buf_offset + k + 2] = smc[i][j].vars[smc_IP3];
		k += grid.num_coupling_species_smc;
	}

	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc * (sendbuf[RIGHT][1] - sendbuf[RIGHT][0] + 1);

	k = 0;
	for (int j = (int) sendbuf[RIGHT][2]; j <= (int) sendbuf[RIGHT][3]; j++) {
		int i = grid.num_ec_circumferentially;
		sendbuf[RIGHT][buf_offset + k + 0] = ec[i][j].vars[ec_Ca];
		sendbuf[RIGHT][buf_offset + k + 1] = ec[i][j].vars[ec_Vm];
		sendbuf[RIGHT][buf_offset + k + 2] = ec[i][j].vars[ec_IP3];
		k += grid.num_coupling_species_ec;
	}
}

/*
 Unpack received data into ghost cells.
 Specific to KNBGR.
*/
void communication_update_recvbuf(grid_parms grid, double** recvbuf, SMC_cell** smc, EC_cell** ec)
{
	int k, buf_offset;

/// UP direction
	buf_offset = grid.added_info_in_send_buf;

	k = 0;
	for (int i = (int) recvbuf[UP][0]; i <= (int) recvbuf[UP][1]; i++) {
		int j = 0;
		smc[i][j].vars[smc_Ca] = recvbuf[UP][buf_offset + k + 0];
		smc[i][j].vars[smc_Vm] = recvbuf[UP][buf_offset + k + 1];
		smc[i][j].vars[smc_IP3] = recvbuf[UP][buf_offset + k + 2];
		k += grid.num_coupling_species_smc;
	}

	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc * (recvbuf[UP][1] - recvbuf[UP][0] + 1);

	k = 0;
	for (int i = (int) recvbuf[UP][2]; i <= (int) recvbuf[UP][3]; i++) {
		int j = 0;
		ec[i][j].vars[ec_Ca] = recvbuf[UP][buf_offset + k + 0];
		ec[i][j].vars[ec_Vm] = recvbuf[UP][buf_offset + k + 1];
		ec[i][j].vars[ec_IP3] = recvbuf[UP][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}

/// DOWN direction
	buf_offset = grid.added_info_in_send_buf;

	k = 0;
	if (grid.flip_array[DOWN] == 0) {
		for (int i = (int) recvbuf[DOWN][0]; i <= (int) recvbuf[DOWN][1];
				i++) {
			int j = grid.num_smc_axially + 1;
			smc[i][j].vars[smc_Ca] = recvbuf[DOWN][buf_offset + k + 0];
			smc[i][j].vars[smc_Vm] = recvbuf[DOWN][buf_offset + k + 1];
			smc[i][j].vars[smc_IP3] = recvbuf[DOWN][buf_offset + k + 2];
			k += grid.num_coupling_species_smc;
		}
	} else if (grid.flip_array[DOWN] == 1) {
		int start = (int) recvbuf[DOWN][0], end = (int) recvbuf[DOWN][1];
		for (int i = end; i >= start; i--) {
			int j = grid.num_smc_axially + 1;
			smc[i][j].vars[smc_Ca] = recvbuf[DOWN][buf_offset + k + 0];
			smc[i][j].vars[smc_Vm] = recvbuf[DOWN][buf_offset + k + 1];
			smc[i][j].vars[smc_IP3] = recvbuf[DOWN][buf_offset + k + 2];
			k += grid.num_coupling_species_smc;
		}
	}

	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc * (recvbuf[DOWN][1] - recvbuf[DOWN][0] + 1);

	k = 0;
	if (grid.flip_array[DOWN] == 0) {
		for (int i = (int) recvbuf[DOWN][2]; i <= (int) recvbuf[DOWN][3];
				i++) {
			int j = grid.num_ec_axially + 1;
			ec[i][j].vars[ec_Ca] = recvbuf[DOWN][buf_offset + k + 0];
			ec[i][j].vars[ec_Vm] = recvbuf[DOWN][buf_offset + k + 1];
			ec[i][j].vars[ec_IP3] = recvbuf[DOWN][buf_offset + k + 2];
			k += grid.num_coupling_species_ec;
		}

	} else if (grid.flip_array[DOWN] == 1) {
		int start = (int) recvbuf[DOWN][2], end = (int) recvbuf[DOWN][3];
		for (int i = end; i >= start; i--) {
			int j = grid.num_ec_axially + 1;
			ec[i][j].vars[ec_Ca] = recvbuf[DOWN][buf_offset + k + 0];
			ec[i][j].vars[ec_Vm] = recvbuf[DOWN][buf_offset + k + 1];
			ec[i][j].vars[ec_IP3] = recvbuf[DOWN][buf_offset + k + 2];
			k += grid.num_coupling_species_ec;
		}
	}

/// LEFT direction
	buf_offset = grid.added_info_in_send_buf;

	k = 0;
	for (int j = (int) recvbuf[LEFT][0]; j <= (int) recvbuf[LEFT][1]; j++) {
		int i = 0;
		smc[i][j].vars[smc_Ca] = recvbuf[LEFT][buf_offset + k + 0];
		smc[i][j].vars[smc_Vm] = recvbuf[LEFT][buf_offset + k + 1];
		smc[i][j].vars[smc_IP3] = recvbuf[LEFT][buf_offset + k + 2];
		k += grid.num_coupling_species_smc;
	}

	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc * (recvbuf[LEFT][1] - recvbuf[LEFT][0] + 1);

	k = 0;
	for (int j = (int) recvbuf[LEFT][2]; j <= (int) recvbuf[LEFT][3]; j++) {
		int i = 0;
		ec[i][j].vars[ec_Ca] = recvbuf[LEFT][buf_offset + k + 0];
		ec[i][j].vars[ec_Vm] = recvbuf[LEFT][buf_offset + k + 1];
		ec[i][j].vars[ec_IP3] = recvbuf[LEFT][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}

/// RIGHT direction
	buf_offset = grid.added_info_in_send_buf;

	k = 0;
	for (int j = (int) recvbuf[RIGHT][0]; j <= (int) recvbuf[RIGHT][1]; j++) {
		int i = grid.num_smc_circumferentially + 1;
		smc[i][j].vars[smc_Ca] = recvbuf[RIGHT][buf_offset + k + 0];
		smc[i][j].vars[smc_Vm] = recvbuf[RIGHT][buf_offset + k + 1];
		smc[i][j].vars[smc_IP3] = recvbuf[RIGHT][buf_offset + k + 2];
		k += grid.num_coupling_species_smc;
	}

	buf_offset = grid.added_info_in_send_buf
			+ grid.num_coupling_species_smc * (recvbuf[RIGHT][1] - recvbuf[RIGHT][0] + 1);

	k = 0;
	for (int j = (int) recvbuf[RIGHT][2]; j <= (int) recvbuf[RIGHT][3]; j++) {
		int i = grid.num_ec_circumferentially + 1;
		ec[i][j].vars[ec_Ca] = recvbuf[RIGHT][buf_offset + k + 0];
		ec[i][j].vars[ec_Vm] = recvbuf[RIGHT][buf_offset + k + 1];
		ec[i][j].vars[ec_IP3] = recvbuf[RIGHT][buf_offset + k + 2];
		k += grid.num_coupling_species_ec;
	}
}
