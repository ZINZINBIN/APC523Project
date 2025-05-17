import argparse, os
from scipy.io import savemat
from src.PIC import PIC
from src.dist import BumpOnTail1D
from src.util import (
    generate_bump_on_tail_figure,
    generate_bump_on_tail_snapshot,
    generate_hamiltonian_analysis,
    generate_bump_on_tail_gif,
    generate_bump_on_tail_dist_gif,
    generate_distribution_snapshot,
    generate_distribution_figure,
    generate_v_distribution_figure,
)

def parsing():
    parser = argparse.ArgumentParser(description="1D Electrostatic Particle-In-Cell code for plasma kinetic simulation")
    parser.add_argument("--num_particle", type = int, default = 40000)
    parser.add_argument("--num_mesh", type = int, default = 1000)
    parser.add_argument("--method", type = str, default = "leapfrog", choices=["midpoint","leapfrog", "verlet", "implicit"])
    parser.add_argument("--solver", type=str, default="Gauss", choices=["SOR", "Gauss"])
    parser.add_argument("--interpol", type = str, default = "TSC", choices=["CIC", "TSC"])
    parser.add_argument("--t_min", type = float, default = 0)
    parser.add_argument("--t_max", type = float, default = 50.0)
    parser.add_argument("--dt", type = float, default = 0.01)
    parser.add_argument("--L", type = float, default = 50)
    parser.add_argument("--n0", type = float, default = 1.0)
    parser.add_argument("--gamma", type = float, default = 5.0)
    parser.add_argument("--A", type=float, default=0.02)
    parser.add_argument("--n_mode", type=int, default=5)
    parser.add_argument("--a", type=float, default=0.2)
    parser.add_argument("--v0", type=float, default=3.0)
    parser.add_argument("--sigma", type=float, default=0.5)
    parser.add_argument("--CFL", type=bool, default=False)
    parser.add_argument("--use_animation", type = bool, default = False)
    parser.add_argument("--plot_freq", type = int, default = 50)
    parser.add_argument("--simcase", type=str, default="bump-on-tail", choices = ["two-stream", "bump-on-tail"])
    args = vars(parser.parse_args())
    return args

if __name__ == "__main__":

    args = parsing()
    
    savedir = os.path.join("./result/bump-on-tail", "{}_{}_N_{}_Nm_{}_dt_{}".format(args["method"], args['interpol'], args['num_particle'], args['num_mesh'], args['dt']))
    
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    # Initial distribution: Bump-On-Tail distribution
    dist = BumpOnTail1D(a = args['a'], v0 = args['v0'], sigma = args['sigma'], n_samples=args['num_particle'], L = args['L'])

    sim = PIC(
        N=args["num_particle"],
        N_mesh=args["num_mesh"],
        method=args["method"],
        solver=args["solver"],
        interpol=args["interpol"],
        n0=args["n0"],
        L=args["L"],
        A=args["A"],
        n_mode=args["n_mode"],
        dt=args["dt"],
        tmin=args["t_min"],
        tmax=args["t_max"],
        gamma=args["gamma"],
        simcase=args["simcase"],
        init_dist=dist,
        CFL = args['CFL']
    )

    snapshot, E, KE, PE = sim.solve()

    h_idx = dist.high_indx

    # plot pic simulation figure
    generate_bump_on_tail_snapshot(snapshot[:,-1], savedir, "{}_snapshot_{}_{}.png".format(args['simcase'], args['interpol'], args['method']), xmin = 0, xmax = args['L'], vmin = -10.0, vmax = 10.0, high_electron_indice=h_idx)
    generate_bump_on_tail_figure(snapshot, savedir, "{}_evolution_{}_{}.png".format(args['simcase'], args['interpol'], args['method']), xmin = 0, xmax = args['L'], vmin = -10.0, vmax = 10.0, high_electron_indice=h_idx)
    generate_hamiltonian_analysis(args['t_max'], E, KE, PE, savedir, "{}_hamiltonian_{}_{}.png".format(args['simcase'], args['interpol'], args['method']))

    # plot distribution
    generate_distribution_snapshot(snapshot[:,-1], savedir, "{}_snapshot_dist_{}_{}.png".format(args['simcase'], args['interpol'], args['method']), xmin = 0, xmax = args['L'], vmin = -10.0, vmax = 10.0)
    generate_distribution_figure(snapshot, savedir, "{}_evolution_dist_{}_{}.png".format(args['simcase'], args['interpol'], args['method']), xmin = 0, xmax = args['L'], vmin = -10.0, vmax = 10.0)
    generate_v_distribution_figure(snapshot, savedir,"{}_evolution_v_dist_{}_{}.png".format(args['simcase'], args['interpol'], args['method']), vmin = -10.0, vmax = 10.0)
    
    if args['use_animation']:
        generate_bump_on_tail_gif(snapshot, savedir, "{}_simulation_{}_{}.gif".format(args['simcase'], args['interpol'], args['method']), 0, args['L'], -10.0, 10.0, args['plot_freq'], h_idx)
        generate_bump_on_tail_dist_gif(snapshot, savedir, "{}_simulation_dist_{}_{}.gif".format(args['simcase'], args['interpol'], args['method']), 0, args['L'], -10.0, 10.0, args['plot_freq'], h_idx)
        
    mdic = {
        "snapshot": snapshot,
        "N": args["num_particle"],
        "N_mesh": args["num_mesh"],
        "n0": args["n0"],
        "L": args["L"],
        "dt": args["dt"],
        "tmin": args["t_min"],
        "tmax": args["t_max"],
        "v0": args['v0'],
        "n_mode":args['n_mode'],
        "A":args['A'],
        "H": E,
    }

    # save data
    if not os.path.exists('./data/bump-on-tail'):
        os.makedirs("./data/bump-on-tail")
    
    savemat(file_name = os.path.join("./data/bump-on-tail", "{}_{}_N_{}_Nm_{}_dt_{}.mat".format(args["method"], args['interpol'], args['num_particle'], args['num_mesh'], args['dt'])), mdict=mdic, do_compression=True)
