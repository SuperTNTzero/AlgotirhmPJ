# 哈希DP算法伪代码

## 1. 哈希函数

```
HASH-FUNCTION(sequence, base_num, mod)
    // 计算序列的哈希值，使用多项式哈希函数
    // 输入：sequence - 待哈希的序列
    //      base_num - 基数（通常为质数，如101）
    //      mod - 模数（通常为大质数，如10^9+7）
    // 输出：序列的哈希值

    BASE-VAL ← 一个映射表，将碱基映射为数值
                {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    hash_val ← 0
    
    for i ← 0 to length[sequence] - 1 do
        base_val ← BASE-VAL[sequence[i]]
        hash_val ← (hash_val * base_num + base_val) mod mod
    
    return hash_val
```

## 2. 滚动哈希

```
ROLLING-HASH(prev_hash, prev_char, next_char, length, base_num, mod)
    // 计算滚动哈希值，用于快速更新子序列的哈希值
    // 输入：prev_hash - 前一个哈希值
    //      prev_char - 要移除的字符（序列最左侧的字符）
    //      next_char - 要添加的字符（序列最右侧的新字符）
    //      length - 子序列长度
    //      base_num - 基数（通常为质数，如101）
    //      mod - 模数（通常为大质数，如10^9+7）
    // 输出：更新后的哈希值

    BASE-VAL ← 一个映射表，将碱基映射为数值
                {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    
    prev_val ← BASE-VAL[prev_char]
    next_val ← BASE-VAL[next_char]
    
    // 计算最高位的贡献值
    highest_digit_val ← (prev_val * POWER(base_num, length-1, mod)) mod mod
    
    // 移除最高位的贡献并添加新字符的贡献
    new_hash ← ((prev_hash - highest_digit_val + mod) * base_num + next_val) mod mod
    
    return new_hash
```

## 3. 预计算幂值

```
PRECOMPUTE-POWERS(base_num, max_length, mod)
    // 预计算幂值，避免重复计算
    // 输入：base_num - 基数
    //      max_length - 最大长度
    //      mod - 模数
    // 输出：幂值数组

    powers ← 一个新数组，大小为max_length
    powers[0] ← 1
    
    for i ← 1 to max_length - 1 do
        powers[i] ← (powers[i-1] * base_num) mod mod
    
    return powers
```

## 4. 处理序列查找重复子序列

```
PROCESS-SEQUENCE(reference, query, is_reversed, min_length, max_length)
    // 处理序列查找重复子序列
    // 输入：reference - 参考序列
    //      query - 查询序列
    //      is_reversed - 是否为反向重复
    //      min_length - 最小子序列长度
    //      max_length - 最大子序列长度
    // 输出：找到的重复序列列表

    results ← 空列表
    unique_sequences ← 空集合
    
    ref_len ← length[reference]
    query_len ← length[query]
    base_num ← 101
    mod ← 10^9 + 7
    
    // 设置最大长度限制
    if max_length = NIL then
        max_length ← min(ref_len, query_len)
    else
        max_length ← min(max_length, ref_len, query_len)
    
    // 预计算幂值
    powers ← PRECOMPUTE-POWERS(base_num, max(ref_len, query_len), mod)
    
    // 对于每个可能的子序列长度
    for length ← min_length to max_length do
        // 忽略单个碱基的重复序列
        if length = 1 then
            continue
        
        // 提前终止条件
        if |results| > 100 and length < min_length + 5 then
            break
        
        // 创建哈希表，存储reference中所有子序列的哈希值
        ref_hash_map ← 空哈希表
        
        // 计算reference中所有子序列的哈希值
        if ref_len ≥ length then
            curr_hash ← HASH-FUNCTION(reference[0..length-1], base_num, mod)
            
            if curr_hash ∉ ref_hash_map then
                ref_hash_map[curr_hash] ← 空列表
            APPEND(ref_hash_map[curr_hash], 0)
            
            // 使用滚动哈希计算其余子序列的哈希值
            for i ← 1 to ref_len - length do
                curr_hash ← ROLLING-HASH(curr_hash, reference[i-1], reference[i+length-1], length, base_num, mod)
                
                if curr_hash ∉ ref_hash_map then
                    ref_hash_map[curr_hash] ← 空列表
                APPEND(ref_hash_map[curr_hash], i)
        
        // 计算query中所有子序列的哈希值并查找重复
        if query_len ≥ length then
            // 创建查询序列的哈希表
            query_hash_map ← 空哈希表
            
            // 计算query中第一个子序列的哈希值
            curr_hash ← HASH-FUNCTION(query[0..length-1], base_num, mod)
            
            if curr_hash ∉ query_hash_map then
                query_hash_map[curr_hash] ← 空列表
            APPEND(query_hash_map[curr_hash], 0)
            
            // 使用滚动哈希计算其余子序列的哈希值
            for i ← 1 to query_len - length do
                curr_hash ← ROLLING-HASH(curr_hash, query[i-1], query[i+length-1], length, base_num, mod)
                
                if curr_hash ∉ query_hash_map then
                    query_hash_map[curr_hash] ← 空列表
                APPEND(query_hash_map[curr_hash], i)
            
            // 对于reference中的每个子序列，检查在query中是否有重复
            for each hash_val in ref_hash_map do
                if hash_val ∈ query_hash_map then
                    // 验证实际序列
                    for each ref_pos in ref_hash_map[hash_val] do
                        ref_subseq ← reference[ref_pos..ref_pos+length-1]
                        
                        // 如果序列已经处理过，则跳过
                        if ref_subseq ∈ unique_sequences then
                            continue
                        
                        query_positions ← 空列表
                        
                        // 查找query中所有匹配的位置
                        for each query_pos in query_hash_map[hash_val] do
                            query_subseq ← query[query_pos..query_pos+length-1]
                            if ref_subseq = query_subseq then
                                APPEND(query_positions, query_pos)
                        
                        // 检查序列在reference中的出现次数
                        ref_all_positions ← 空列表
                        ref_check_pos ← 0
                        while true do
                            ref_check_pos ← FIND(reference, ref_subseq, ref_check_pos)
                            if ref_check_pos = -1 then
                                break
                            APPEND(ref_all_positions, ref_check_pos)
                            ref_check_pos ← ref_check_pos + 1
                        
                        // 只有当序列在reference中仅出现一次，但在query中出现多次时，才认为是重复序列
                        if |ref_all_positions| = 1 and |query_positions| > 1 then
                            // 计算额外重复次数
                            repeat_count ← |query_positions| - 1
                            
                            // 添加重复结果
                            ADD-REPEAT-RESULT(results, ref_subseq, [ref_pos], query_positions, is_reversed, repeat_count, length)
                            
                            // 将序列添加到已处理集合中
                            ADD(unique_sequences, ref_subseq)
    
    return results
```

## 5. 添加重复变异结果

```
ADD-REPEAT-RESULT(results, subseq, ref_positions, query_positions, is_reversed, repeat_count, length)
    // 添加重复变异结果
    // 输入：results - 结果列表
    //      subseq - 重复序列
    //      ref_positions - 在reference中的位置列表
    //      query_positions - 在query中的位置列表
    //      is_reversed - 是否为反向重复
    //      repeat_count - 重复次数（额外出现的次数）
    //      length - 序列长度
    // 输出：无，直接修改results列表

    // 检查是否已经有相同序列的结果
    for each result in results do
        if result.sequence = subseq then
            // 序列已存在，不再添加
            return
    
    // 检查是否已经有更长的序列包含了这个子序列
    is_contained ← false
    for each result in results do
        if result.reversed = is_reversed and
           length[result.sequence] > length and
           ∀p ∈ ref_positions: p ∈ result.ref_positions then
            is_contained ← true
            break
    
    if not is_contained then
        new_result ← 创建新的结果对象
        new_result.sequence ← subseq
        new_result.ref_positions ← ref_positions
        new_result.query_positions ← query_positions
        new_result.repeat_count ← repeat_count
        new_result.reversed ← is_reversed
        
        APPEND(results, new_result)
```

## 6. 主算法：查找重复序列

```
FIND-REPEATS-HASH-DP(reference, query, min_length, max_length, use_parallel)
    // 使用哈希和动态规划方法查找重复序列及其位置
    // 输入：reference - 参考序列
    //      query - 查询序列
    //      min_length - 最小子序列长度
    //      max_length - 最大子序列长度
    //      use_parallel - 是否使用并行处理
    // 输出：找到的重复序列列表

    unique_sequences ← 空集合
    results ← 空列表
    
    // 设置最大长度限制
    if max_length = NIL then
        max_length ← min(length[reference], length[query])
    else
        max_length ← min(max_length, length[reference], length[query])
    
    // 如果序列太短，直接返回空结果
    if min_length > max_length then
        return results
    
    // 使用并行处理同时处理正向和反向序列
    if use_parallel and max_length > 20 then
        // 并行处理正向序列和反向互补序列
        forward_results ← 并行执行 PROCESS-SEQUENCE(reference, query, false, min_length, max_length)
        
        query_rev ← REVERSE-COMPLEMENT(query)
        reverse_results ← 并行执行 PROCESS-SEQUENCE(reference, query_rev, true, min_length, max_length)
        
        // 合并结果
        APPEND(results, forward_results)
        APPEND(results, reverse_results)
    else
        // 处理正向序列
        PROCESS-SEQUENCE(reference, query, false, results, min_length, max_length, unique_sequences)
        
        // 处理反向互补序列
        query_rev ← REVERSE-COMPLEMENT(query)
        PROCESS-SEQUENCE(reference, query_rev, true, results, min_length, max_length, unique_sequences)
    
    // 按序列长度降序排序
    SORT(results, 按照序列长度降序)
    
    return results
```

## 7. 反向互补序列生成

```
REVERSE-COMPLEMENT(sequence)
    // 生成DNA序列的反向互补序列
    // 输入：sequence - DNA序列
    // 输出：反向互补序列

    complement ← {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    result ← 空字符串
    
    for i ← length[sequence] - 1 downto 0 do
        base ← sequence[i]
        if base ∈ complement then
            result ← result + complement[base]
        else
            result ← result + base
    
    return result
```

## 算法复杂度分析

- **时间复杂度**：
  - 哈希函数计算：O(L)，其中L是序列长度
  - 滚动哈希更新：O(1)
  - 处理序列查找重复：O(N * M * L)，其中N是参考序列长度，M是查询序列长度，L是最大子序列长度
  - 使用并行处理可以将时间复杂度降低为O(N * M * L / P)，其中P是处理器数量

- **空间复杂度**：
  - 哈希表存储：O(N + M)，其中N是参考序列长度，M是查询序列长度
  - 结果存储：O(K)，其中K是找到的重复序列数量

## 算法优化策略

1. **滚动哈希**：使用滚动哈希算法，避免重复计算子序列哈希值
2. **预计算幂值**：预先计算所有可能用到的幂值，避免重复计算
3. **并行处理**：同时处理正向和反向互补序列，提高效率
4. **提前终止**：当找到足够多的长序列时，可以跳过短序列的处理
5. **哈希冲突处理**：使用字符串比较验证哈希匹配，避免哈希冲突导致的错误结果

## 运行结果

![运行结果](D:/AlgorithmPJ/HashAlgoResult.png)

## 时间复杂度

O(n * n)

## Github仓库：https://github.com/SuperTNTzero/AlgotirhmPJ/tree/master
