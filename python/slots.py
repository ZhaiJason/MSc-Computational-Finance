def jump(nums) -> int:

        end_idx = len(nums) - 1
        pre_idx = end_idx
        jumps = 0

        while pre_idx > 0:
            for idx in range(end_idx - 1, -1, -1):
                if idx + nums[idx] >= end_idx:
                    pre_idx = idx
            jumps += 1
            end_idx = pre_idx
        
        return jumps

print(jump([2,3,1,1,4]))