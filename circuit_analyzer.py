"""
Circuit Current Analyzer
Solves for currents through resistors in an electric circuit using nodal analysis.

Input file format (circuit.txt):
- Resistors: two letters = resistance in ohms (e.g., "ab = 100")
- Voltage source node: single letter = voltage (e.g., "a = 5")
- Ground/reference node: single letter = 0 (e.g., "d = 0")

Output: currents through each resistor in amperes
"""

import numpy as np
import re
from collections import defaultdict


def parse_circuit_file(filename):
    """
    Parse the circuit.txt file and extract resistors and node voltages.
    
    Returns:
        resistors: dict of {(node1, node2): resistance}
        fixed_voltages: dict of {node: voltage}
    """
    resistors = {}
    fixed_voltages = {}
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Remove comments and extra whitespace
    lines = content.strip().split('\n')
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        # Parse line: "xy = value" format
        match = re.match(r'^([a-zA-Z]+)\s*=\s*([\d.]+)$', line)
        if match:
            identifier = match.group(1).lower()
            value = float(match.group(2))
            
            if len(identifier) == 2:
                # It's a resistor between two nodes
                node1, node2 = identifier[0], identifier[1]
                # Store with sorted tuple to avoid duplicates
                key = tuple(sorted([node1, node2]))
                resistors[key] = value
            elif len(identifier) == 1:
                # It's a fixed voltage at a node
                fixed_voltages[identifier] = value
            else:
                print(f"Warning: Ignoring invalid identifier '{identifier}'")
    
    return resistors, fixed_voltages


def get_all_nodes(resistors, fixed_voltages):
    """Extract all unique nodes from the circuit."""
    nodes = set()
    for (n1, n2) in resistors.keys():
        nodes.add(n1)
        nodes.add(n2)
    for node in fixed_voltages.keys():
        nodes.add(node)
    return sorted(list(nodes))


def solve_circuit(resistors, fixed_voltages):
    """
    Solve the circuit using nodal analysis (Modified Nodal Analysis).
    
    Uses Kirchhoff's Current Law: sum of currents at each node = 0
    And Ohm's Law: I = V / R
    
    Returns:
        node_voltages: dict of {node: voltage}
        currents: dict of {(node1, node2): current} (positive = flow from node1 to node2)
    """
    nodes = get_all_nodes(resistors, fixed_voltages)
    n = len(nodes)
    
    if n == 0:
        return {}, {}
    
    # Create node index mapping
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    
    # Check for ground node
    ground_node = None
    for node, voltage in fixed_voltages.items():
        if voltage == 0:
            ground_node = node
            break
    
    if ground_node is None:
        # If no ground specified, use the first fixed voltage node as reference
        if fixed_voltages:
            ground_node = list(fixed_voltages.keys())[0]
        else:
            raise ValueError("Circuit must have at least one node with fixed voltage (ground)")
    
    # Build conductance matrix G and current vector I
    # Using nodal analysis: G * V = I
    # For each node: sum of (V_node - V_neighbor) / R = I_injected
    
    G = np.zeros((n, n))
    I = np.zeros(n)
    
    # Build conductance matrix from resistors
    for (n1, n2), resistance in resistors.items():
        if resistance <= 0:
            raise ValueError(f"Invalid resistance {resistance} between {n1} and {n2}")
        
        conductance = 1.0 / resistance
        i, j = node_to_idx[n1], node_to_idx[n2]
        
        # Self conductance (diagonal)
        G[i, i] += conductance
        G[j, j] += conductance
        
        # Mutual conductance (off-diagonal)
        G[i, j] -= conductance
        G[j, i] -= conductance
    
    # Apply fixed voltage constraints
    # For nodes with fixed voltages, replace their equation with V_node = V_fixed
    for node, voltage in fixed_voltages.items():
        idx = node_to_idx[node]
        # Clear the row
        G[idx, :] = 0
        # Set diagonal to 1
        G[idx, idx] = 1
        # Set the voltage
        I[idx] = voltage
    
    # Solve the system
    try:
        V = np.linalg.solve(G, I)
    except np.linalg.LinAlgError:
        raise ValueError("Circuit has no unique solution. Check connections and constraints.")
    
    # Create voltage dictionary
    node_voltages = {node: V[node_to_idx[node]] for node in nodes}
    
    # Calculate currents through each resistor using Ohm's law
    # Current flows from higher voltage to lower voltage
    currents = {}
    for (n1, n2), resistance in resistors.items():
        v1 = node_voltages[n1]
        v2 = node_voltages[n2]
        # Current from n1 to n2 (positive if v1 > v2)
        current = (v1 - v2) / resistance
        currents[(n1, n2)] = current
    
    return node_voltages, currents


def write_output(output_filename, currents, node_voltages=None):
    """Write the calculated currents to an output file."""
    with open(output_filename, 'w') as f:
        f.write("# Circuit Analysis Results\n")
        f.write("# Currents through resistors (in Amperes)\n")
        f.write("# Positive value: current flows from first node to second node\n\n")
        
        if node_voltages:
            f.write("# Node Voltages:\n")
            for node in sorted(node_voltages.keys()):
                f.write(f"# {node} = {node_voltages[node]:.6f} V\n")
            f.write("\n")
        
        f.write("# Resistor Currents:\n")
        for (n1, n2), current in sorted(currents.items()):
            f.write(f"{n1}{n2} = {current:.6f}\n")


def main():
    """Main function to run the circuit analyzer."""
    input_file = "circuit.txt"
    output_file = "output.txt"
    
    print(f"Reading circuit from '{input_file}'...")
    
    try:
        resistors, fixed_voltages = parse_circuit_file(input_file)
    except FileNotFoundError:
        print(f"Error: Could not find '{input_file}'")
        print("Please create a circuit.txt file with your circuit definition.")
        return
    
    print(f"Found {len(resistors)} resistor(s) and {len(fixed_voltages)} fixed voltage node(s)")
    
    if not resistors:
        print("Error: No resistors found in the circuit.")
        return
    
    if len(fixed_voltages) < 2:
        print("Error: Need at least 2 fixed voltage nodes (source and ground).")
        return
    
    # Check for ground
    has_ground = any(v == 0 for v in fixed_voltages.values())
    if not has_ground:
        print("Warning: No ground node (voltage = 0) specified.")
    
    print("\nSolving circuit using nodal analysis...")
    
    try:
        node_voltages, currents = solve_circuit(resistors, fixed_voltages)
    except ValueError as e:
        print(f"Error solving circuit: {e}")
        return
    
    print("\nNode Voltages:")
    for node in sorted(node_voltages.keys()):
        print(f"  {node} = {node_voltages[node]:.6f} V")
    
    print("\nResistor Currents:")
    for (n1, n2), current in sorted(currents.items()):
        abs_current = abs(current)
        direction = f"{n1}->{n2}" if current >= 0 else f"{n2}->{n1}"
        print(f"  {n1}{n2} = {current:.6f} A (flows {direction})")
    
    write_output(output_file, currents, node_voltages)
    print(f"\nResults written to '{output_file}'")


if __name__ == "__main__":
    main()
