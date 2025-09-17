#!/usr/bin/env python3
"""
Cost estimation script for Donut Falls pipeline on AWS
Provides rough estimates for different configurations
"""

import argparse
import sys

def estimate_costs(num_samples, profile='standard', instance_type='standard'):
    """
    Estimate costs based on sample count and configuration
    
    Rough estimates based on AWS pricing (us-east-1)
    Actual costs will vary based on region, instance types, and usage patterns
    """
    
    # Base cost per sample estimates (USD) - AWS pricing
    cost_estimates = {
        'standard': {
            'flye': 3.50,
            'unicycler': 2.80,
            'raven': 2.20
        },
        'cost_optimized': {
            'flye': 1.75,
            'unicycler': 1.40,
            'raven': 1.10
        },
        'aws_optimized': {
            'flye': 1.20,
            'unicycler': 0.95,
            'raven': 0.80
        },
        'unicycler_optimized': {
            'unicycler': 0.40
        }
    }
    
    # Instance type multipliers
    instance_multipliers = {
        'standard': 1.0,
        'spot': 0.3,  # ~70% discount
        'compute_optimized': 0.9,
        'memory_optimized': 1.2
    }
    
    if profile not in cost_estimates:
        print(f"Error: Unknown profile '{profile}'")
        sys.exit(1)
    
    if instance_type not in instance_multipliers:
        print(f"Error: Unknown instance type '{instance_type}'")
        sys.exit(1)
    
    print(f"\n=== Cost Estimation for Donut Falls on AWS ===")
    print(f"Samples: {num_samples}")
    print(f"Profile: {profile}")
    print(f"Instance Type: {instance_type}")
    print(f"{'='*45}")
    
    total_cost = 0
    
    for assembler, cost_per_sample in cost_estimates[profile].items():
        adjusted_cost = cost_per_sample * instance_multipliers[instance_type]
        assembler_total = adjusted_cost * num_samples
        total_cost += assembler_total
        
        print(f"{assembler:15} ${adjusted_cost:6.2f}/sample × {num_samples:2d} = ${assembler_total:7.2f}")
    
    print(f"{'='*45}")
    print(f"{'Total Estimated':15} ${total_cost:7.2f}")
    
    # Show potential savings
    if profile != 'standard':
        standard_cost = cost_estimates['standard']['unicycler'] * instance_multipliers[instance_type] * num_samples
        savings = standard_cost - total_cost
        savings_percent = (savings / standard_cost) * 100
        print(f"{'Potential Savings':15} ${savings:7.2f} ({savings_percent:.1f}%)")
    
    print(f"\nNote: These are rough estimates based on AWS us-east-1 pricing.")
    print("Actual costs depend on:")
    print("- Exact instance types and regions used")
    print("- Data transfer and storage costs")
    print("- Actual runtime and resource utilization")
    print("- AWS pricing changes and spot price fluctuations")
    
    # Show recommended command
    if profile == 'unicycler_optimized':
        print(f"\nRecommended command:")
        print(f"nextflow run UPHL-BioNGS/Donut_Falls \\")
        print(f"  -profile singularity,unicycler_optimized \\")
        print(f"  --sample_sheet sample_sheet.csv \\")
        print(f"  --assembler unicycler")

def main():
    parser = argparse.ArgumentParser(description='Estimate costs for Donut Falls pipeline on AWS')
    parser.add_argument('samples', type=int, help='Number of samples to process')
    parser.add_argument('--profile', choices=['standard', 'cost_optimized', 'aws_optimized', 'unicycler_optimized'], 
                       default='standard', help='Configuration profile to use')
    parser.add_argument('--instance-type', choices=['standard', 'spot', 'compute_optimized', 'memory_optimized'],
                       default='standard', help='Instance type category')
    
    args = parser.parse_args()
    
    estimate_costs(args.samples, args.profile, args.instance_type)

if __name__ == '__main__':
    main()