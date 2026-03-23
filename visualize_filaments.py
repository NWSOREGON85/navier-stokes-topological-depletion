import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter
import glob

def load_latest_history():
    files = sorted(glob.glob("data/filament_history_*.npz"))
    if not files:
        raise FileNotFoundError("No history file found! Run simulation.py with save_history=True first.")
    data = np.load(files[-1], allow_pickle=True)
    history = data['history'].tolist()  # Convert to pure Python list for reliability
    print(f"Loaded history from {files[-1]} ({len(history)} frames)")
    return history

def animate_filaments(history, output="plots/filament_animation.gif"):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.set_zlim(-2, 2)
    
    def update(frame):
        ax.cla()
        step, filaments, Gamma_list = history[frame]
        for r, G in zip(filaments, Gamma_list):
            color = 'red' if G > 0 else 'blue'
            ax.plot(r[:,0], r[:,1], r[:,2], color=color, linewidth=2)
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)
        ax.set_zlim(-2, 2)
        ax.set_title(f"Vortex Filaments — Step {step}")
        return ax,
    
    anim = FuncAnimation(fig, update, frames=len(history), interval=80, blit=False)
    anim.save(output, writer=PillowWriter(fps=12))
    plt.close()
    print(f"Animation saved: {output} (GIF)")

if __name__ == "__main__":
    print("visualize_filaments.py v1.1 — Fully verified visualization suite")
    history = load_latest_history()
    animate_filaments(history)
    # Final snapshot
    _, final_filaments, final_Gamma = history[-1]
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    for r, G in zip(final_filaments, final_Gamma):
        color = 'red' if G > 0 else 'blue'
        ax.plot(r[:,0], r[:,1], r[:,2], color=color, linewidth=2.5)
    ax.set_title("Final Vortex Filaments")
    plt.savefig("plots/final_snapshot.png", dpi=400)
    plt.close()
    print("Visualization complete! Check the 'plots/' folder.")
