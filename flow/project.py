"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
import os
import numpy as np
import mbuild as mb


class CDS(FlowProject):
    pass


class Borah(DefaultSlurmEnvironment):
    hostname_pattern = "borah"
    template = "borah.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="shortgpu",
            help="Specify the partition to submit to."
        )
        parser.add_argument(
            "--exclude",
            default="gpu105",
            help="Specify the type of nodes to exclude."
        )


class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="batch",
            help="Specify the partition to submit to."
        )


@CDS.label
def initial_run_done(job):
    return job.doc.runs >= 1


@CDS.label
def equilibrated(job):
    return job.doc.equilibrated


def make_system(job):
    class CopperLattice(mb.Compound):
        def __init__(self, x, y, z, pillar_height_percent=1/3):
            '''Make sure x and y values are divisible '''
            self.php = pillar_height_percent
            super(CopperLattice, self).__init__()
            spacings = [0.36149, 0.36149, 0.36149]
            angles = [90, 90, 90]
            points = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
            fcc_lattice = mb.Lattice(lattice_spacing=spacings, angles=angles, lattice_points={'A' : points})
            # define Compound that is a copper atom
            cu = mb.Compound(name='Cu')
            # populate lattice with compounds
            self.add(fcc_lattice.populate(compound_dict={'A' : cu}, x=x, y=y, z=z))

        @property
        def Lx(self):
            return self.box.Lx

        @property
        def Ly(self):
            return self.box.Ly

        @property
        def Lz(self):
            return self.box.Lz

        @property
        def slab_z_chunk(self):
            return (1-self.php) * self.Lz

        @property
        def surface_thickness(self):
            pass

        def remove_region(self, x_range, y_range, z_range):
            remove_particles = []
            for p in self.particles():
                if p.pos[0] >= x_range[0] and p.pos[0] < x_range[1]:
                    if p.pos[2] >= z_range[0] and p.pos[2] < z_range[1]:
                        remove_particles.append(p)
            for p in remove_particles:
                self.remove(p)

        def make_two_pillars(self):
            """Remove 3 chunks, leaving 2 pillars"""
            # Left-most region
            x_range = (0, self.Lx / 5)
            x_range2 = (2*self.Lx/5, 3*self.Lx/5)
            x_range3 = (4*self.Lx/5, self.Lx)
            y_range = (0, self.Ly)
            z_range = (self.slab_z_chunk/2, self.Lz - (self.slab_z_chunk/2))
            self.remove_region(x_range, y_range, z_range)
            self.remove_region(x_range2, y_range, z_range)
            self.remove_region(x_range3, y_range, z_range)

        def make_four_pillars(self):
            """Remove 5 chunks, leaving 4 pillars"""
            # Left-most region
            x_range = (0, self.Lx / 9)
            x_range2 = (2*self.Lx/9, 3*self.Lx/9)
            x_range3 = (4*self.Lx/9, 5*self.Lx/9)
            x_range4 = (6*self.Lx/9, 7*self.Lx/9)
            x_range5 = (8*self.Lx/9, self.Lx)
            y_range = (0, self.Ly)
            z_range = (self.slab_z_chunk/2, self.Lz - (self.slab_z_chunk/2))
            self.remove_region(x_range, y_range, z_range)
            self.remove_region(x_range2, y_range, z_range)
            self.remove_region(x_range3, y_range, z_range)
            self.remove_region(x_range4, y_range, z_range)
            self.remove_region(x_range5, y_range, z_range)

    system = CopperLattice(x=job.sp.x, y=job.sp.y, z=job.sp.z, pillar_height_percent=job.sp.php)
    if job.sp.n_pillars == 2:
        system.make_two_pillars()
    if job.sp.n_pillars == 4:
        system.make_four_pillars()
    system.save(job.fn('init.gsd'))

@CDS.post(initial_run_done)
@CDS.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="first-run"
)
def initial_run(job):
    import hoomd
    from hoomd import md
    from hoomd import metal 
	
	make_system(job)
	hoomd_args = f"--single-mpi --mode=cpu"
	hoomd.context.initialize(hoomd_args)

	# Create a 10x10x10 simple cubic lattice of particles with type name A
	hoomd.init.read_gsd(job.fn('init.gsd'))

	nl = md.nlist.cell()
	eam = metal.pair.eam(file=job.project.fn('cu.eam.alloy'), type='Alloy', nlist=nl)
	nl.r_cut = job.sp.r_cut

	# Set up GSD writer:
	hoomd.dump.gsd(
		filename="trajectory.gsd",
		period=job.sp.gsd_period,
		group=hoomd.group.all(),
		overwrite=True,
		phase=0,
		dynamic=["momentum"],
	)
	md.integrate.mode_standard(dt=job.sp.dt)
	# Integrate at constant temperature
	hoomd.md.integrate.nvt(group=hoomd.group.all(), kT=job.sp.kT, tau=job.sp.dt * job.sp.tau)
	# Run for 10,000 time steps
	hoomd.run(job.sp.n_steps)






@CDS.pre(initial_run_done)
@CDS.post(equilibrated)
@CDS.operation(
        directives={"ngpu": 1, "executable": "python -u"},
        name="run-longer"
)
def run_longer(job):
	pass



if __name__ == "__main__":
    CDS(environment=Fry).main()
